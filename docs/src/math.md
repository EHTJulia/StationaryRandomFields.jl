# Mathematical conventions

## Latent prior

For `M = prod(size(latent))` independent elements:

```math
\log p_Z(z)=\sum_{i=1}^M\log p(z_i).
```

For a Laplace marginal with location `μ` and scale `b`:

```math
\log p_Z(z)=-M\log(2b)-\frac{1}{b}\sum_i|z_i-\mu|.
```

The scale is not a standard deviation: its variance is `2b^2`. At the Laplace cusp,
AD uses the zero subgradient. Finite-difference checks avoid the cusp. An image-space
Laplace prior penalizes image values; a wavelet-coefficient prior penalizes coefficients.
They express different models.

## Hartley synthesis and discrete power

For zero-based grid coordinates, define the true multidimensional Hartley matrix:

```math
H_{nk}=\operatorname{cas}\left(2\pi\sum_d\frac{n_d k_d}{m_d}\right),
\qquad Q=H/\sqrt{M},\qquad Q^\mathsf{T}Q=I.
```

This is one Hartley kernel of the summed phase, not a product of one-dimensional
Hartley kernels. `HartleyTransform` applies `Q`. The companion package's unnormalized
`fht` is divided by `sqrt(M)`; using `ifht` without compensating normalization would
produce a different variance.

For centered, unit-variance Gaussian latent coefficients:

```math
x=\mu_x+Q\,\operatorname{diag}(\sqrt{P_k})z,\qquad
C=Q\,\operatorname{diag}(P_k)Q^\mathsf{T}.
```

For latent variance `v`, the covariance is multiplied by `v`. Standard-normal latent
coefficients are recommended to avoid degeneracy between latent scale and spectral
amplitude. A nonzero latent mean is transformed by the same operator; it is not a
constant signal mean. Add a constant signal mean explicitly in the observation model.

`P` is discrete power, or covariance eigenvalue, rather than a continuous spectral
density. With an unnormalized DFT, the centered periodogram convention is
`abs2(fft(x)) / M`. There is no additional factor of 1/2. Converting a continuous density
to discrete power requires frequency-cell volumes and the chosen Fourier convention;
this conversion is not silently inferred.

Frequencies use cycles per coordinate unit: `k / (m_d * spacing_d)`, in unshifted FFT
ordering. Even-size Nyquist bins use the negative representative. On a finite grid,
`SpectralScaling` averages `P(k)` and `P(-k)` using discrete conjugate indices before
computing an amplitude. This also handles Nyquist aliases after rotation/stretching.
`power_map` remains a direct pointwise map; the paired effective power belongs to the
spectral transform.

`SinglePowerLaw` has an explicit DC power, defaulting to zero. `SaturatedPowerLaw` is
flat below a positive inner frequency and zero at/above an outer frequency. Zero-power
modes require no division. Their amplitude derivative is defined to be zero at exactly
zero power; changing the active band at a hard cutoff remains nonsmooth.

## Spectrum composition

`p + q` adds powers. Scaling power by `c >= 0` scales signal amplitude by `sqrt(c)`.
Scaling a signal by `a` scales its power by `abs2(a)`. Summing signals made from the same
latent sample introduces cross terms and does not implement `p + q`.

`stretched(p, s...)` evaluates `p(ν ./ s)`. It changes frequency-coordinate widths and
does not normalize integrated power. `rotated(p, R)` evaluates `p(R * ν)` for a proper
orthogonal matrix. `rotated(p, angle)` uses the corresponding two-dimensional matrix.
Modifiers compose in application order; the last coordinate modifier is evaluated first
when pulling a frequency back to the base model.

## Wavelets

The initial implementation uses periodic orthogonal filter-bank DWT packing compatible
with Wavelets.jl. The approximation occupies the low-index corner; each subsequent
level transforms that corner. Every dimension must be divisible by `2^levels`.
CPU ranks 1–3 use Wavelets.jl. A separable implementation extends the same convention to
higher ranks, CUDA arrays, and Reactant. Float32 filters are converted once at preparation.

For synthesis `x = S D z`, coefficient weights belong to `D`; a Laplace latent prior
still evaluates `z`. For orthogonal `S` and invertible `D`, this corresponds, up to constants
for fixed weights, to an L1 penalty on `D^{-1} S^T x`. Redundant/custom synthesis does not
in general define that analysis penalty on `x`.

## Inference

```math
\log p(z,\theta\mid y)=\log p(y\mid T_\theta(z))+
\log p_Z(z)+\log p(\theta)+C.
```

The sampled variables are `z` and model parameters. No signal density or Jacobian term
is evaluated. `logpdf` on a transformed field raises an informative error.
