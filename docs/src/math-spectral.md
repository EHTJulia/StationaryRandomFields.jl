# Spectral transforms

Gaussian latent coefficients and spectral scaling produce a Gaussian signal
with the chosen discrete power spectrum.

## Hartley synthesis and discrete power

For a grid of shape ``(m_1,\ldots,m_N)`` with ``M=\prod_d m_d`` elements
and zero-based indices, define the multidimensional Hartley matrix:

```math
H_{nk}=\operatorname{cas}\left(2\pi\sum_d\frac{n_d k_d}{m_d}\right),
\qquad Q=H/\sqrt{M},\qquad Q^\mathsf{T}Q=I.
```

This is one Hartley kernel of the summed phase, not a product of one-dimensional
Hartley kernels. `HartleyTransform` applies ``Q``.

For centered, unit-variance Gaussian latent coefficients:

```math
x=\mu_x+Q\,\operatorname{diag}(\sqrt{P_k})z,\qquad
C=Q\,\operatorname{diag}(P_k)Q^\mathsf{T}.
```

For latent variance ``v``, the covariance is multiplied by ``v``. Standard-normal latent
coefficients are recommended to avoid degeneracy between latent scale and spectral
amplitude. A nonzero latent mean is transformed by the same operator; it is not a
constant signal mean. Add a constant signal mean explicitly in the observation model.

``P`` is discrete power, or covariance eigenvalue, rather than a continuous spectral
density. With an unnormalized DFT, the centered periodogram convention is
``|\operatorname{DFT}(x)_k|^2/M``. Converting a continuous density
to discrete power requires frequency-cell volumes and the chosen Fourier convention.

Frequencies use cycles per coordinate unit: ``f_{k,d}=k_d/(m_d\Delta_d)``, in unshifted FFT
ordering. Even-size Nyquist bins use the negative representative. On a finite grid,
`SpectralScaling` averages ``P(k)`` and ``P(\bar k)`` using discrete conjugate indices before
computing an amplitude. This also handles Nyquist aliases after rotation/stretching.
`power_map` remains a direct pointwise map; the paired effective power belongs to the
spectral transform.

## Spectrum composition

`p + q` adds powers. Scaling power by ``c\geq 0`` scales signal amplitude by ``\sqrt{c}``.
Scaling a signal by ``a`` scales its power by ``|a|^2``. Summing signals made from the same
latent sample introduces cross terms and does not implement `p + q`.

`stretched(p, s...)` evaluates ``p(\nu_1/s_1,\ldots,\nu_N/s_N)``. It changes frequency-coordinate widths and
does not normalize integrated power. `rotated(p, R)` evaluates ``p(R\nu)`` for a proper
orthogonal matrix. `rotated(p, angle)` uses the corresponding two-dimensional matrix.
Modifiers compose in application order; the last coordinate modifier is evaluated first
when pulling a frequency back to the base model.

## Predefined spectral shapes

Let ``r=\lVert f\rVert`` be radial frequency in cycles per coordinate unit,
``k=2\pi r``, and ``N`` the number of dimensions.

- **`ConstantPowerSpectrum{N}(a)`** — flat power:

  ```math
  P(r)=a,\qquad a\geq 0.
  ```

- **`SinglePowerLaw{N}(beta; dc=0)`** — a power law with a separately specified DC mode:

  ```math
  P(r)=\begin{cases}r^\beta,&r>0,\\ \mathrm{dc},&r=0.\end{cases}
  ```

- **`SaturatedPowerLaw{N}(beta, inscale, outscale)`** — flat below ``f_{\rm in}``,
  power-law decay for negative ``\beta``, and a hard outer cutoff:

  ```math
  P(r)=\begin{cases}
  \max(r,f_{\rm in})^\beta,&r<f_{\rm out},\\
  0,&r\geq f_{\rm out},
  \end{cases}
  \qquad 0<f_{\rm in}<f_{\rm out}.
  ```

- **`MaternPS{N}(rho, nu)`** — positive range ``\rho`` and smoothness ``\nu``:

  ```math
  P(r)=\left(1+\frac{k^2}{\kappa^2}\right)^{-\nu-N/2},
  \qquad \kappa=\frac{\sqrt{8\nu}}{\rho}.
  ```

- **`MarkovPS{N}((rho1, rho2, ...))`** — inverse-polynomial power;
  one scale gives a Lorentzian shape:

  ```math
  P(r)=\left[1+\sum_{j=1}^{J}(\rho_j k)^{2j}\right]^{-1}.
  ```

  Scales are nonnegative, with at least one positive entry. The tuple length sets
  the polynomial order.

Matérn and Markov have unit DC power, not unit marginal variance. Use `renormed`
to change the power scale. These are finite-grid spectral models; the Markov
shape does not by itself imply a sparse lattice precision matrix.

### Shape comparison

The curves below use two-dimensional models and are scaled to agree at
``r_0=0.04`` so their shapes can be compared independently of amplitude.

```@example spectral_shapes
using StationaryRandomFields, CairoMakie

models = (
    ConstantPowerSpectrum{2}(1.0),
    SinglePowerLaw{2}(-3.0),
    SaturatedPowerLaw{2}(-3.0, 0.035, 0.18),
    MaternPS{2}(8.0, 1.5),
    MarkovPS{2}((8.0,)),
)
labels = (
    "Constant", "Single power law: beta = -3", "Saturated: 0.035–0.18",
    "Matérn: rho = 8, nu = 1.5", "Markov: rho = 8",
)
frequency = 10 .^ range(-3, log10(0.5); length = 600)
fig = Figure(; size = (950, 520), fontsize = 16)
ax = Axis(
    fig[1, 1]; xlabel = "Radial frequency [cycles / unit]",
    ylabel = "Relative power P(r) / P(0.04)", xscale = log10, yscale = log10
)
for (model, label) in zip(models, labels)
    local power = [power_point(model, f, 0.0) / power_point(model, 0.04, 0.0) for f in frequency]
    # Zero power has no position on a logarithmic axis.
    lines!(ax, frequency, map(p -> p > 0 ? p : NaN, power); label, linewidth = 2)
end
vlines!(ax, [0.18]; color = :gray, linestyle = :dash)
axislegend(ax; position = :lb)
save("spectral-shapes.svg", fig)
nothing # hide
```

![Comparison of the five predefined spectral shapes](spectral-shapes.svg)

The dashed line marks the saturated model's outer cutoff: its power is exactly
zero beyond that line. DC is omitted from the logarithmic plot.
For their image-space textures, see [Gaussian image textures](tutorial-gaussian-2d.md).

## Periodic covariance

The sampler and covariance use the same `FieldGrid` and nonnegative power model.
This API describes a **finite periodic grid**, rather than a continuous random field
specified by an integral over all frequencies.

For grid shape ``(n_1,\ldots,n_N)``, spacing ``\Delta_d``, and ``M=\prod_d n_d``,
the periods are ``L_d=n_d\Delta_d``. Frequencies are measured in cycles per physical
coordinate unit and use the FFT representatives documented in [Mathematical Basis](math.md).
Let ``\bar k`` be the lattice index for the conjugate mode, including Nyquist aliases.
The effective power is

```math
\lambda_k=\frac{P(\nu_k)+P(\nu_{\bar k})}{2}.
```

For independent Gaussian latent coefficients with variance ``v``, the covariance is

```math
K(x,y)=\frac{v}{M}\sum_k\lambda_k
\cos\bigl(2\pi\nu_k\cdot(x-y)\bigr).
```

The factor ``1/M`` follows from the orthonormal Hartley synthesis. There is no
frequency-bin volume factor or continuous-spectrum integration. The latent location
affects the signal mean and is excluded from this covariance.

At grid positions ``x_d=j_d\Delta_d``, this kernel exactly matches the Gaussian/FHT
sampler covariance, up to floating-point rounding. The implementation averages the
same **lattice** conjugate pair as `SpectralScaling`; substituting `P(-ν)` can give a
different result on even grids with anisotropic spectra and Nyquist aliases.

For arbitrary real coordinates, the formula defines a finite periodic cosine kernel.
It has a nonnegative sine/cosine feature representation and is positive semidefinite.
This specifies an off-grid extension of the grid covariance; it is not a claim that
the original Hartley latent coefficients alone define those off-grid samples.
In particular, the chosen Nyquist representatives specify this extension and the grid
covariance does not uniquely determine an off-grid model.
