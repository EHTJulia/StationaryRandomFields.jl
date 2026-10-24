# Finite periodic covariance kernels

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

## KernelFunctions extension

Load KernelFunctions to activate the optional package extension. The wrapper follows
the [official custom-kernel interface](https://juliagaussianprocesses.github.io/KernelFunctions.jl/stable/create_kernel/):
it subtypes `KernelFunctions.Kernel` and implements scalar pair evaluation.

```julia
using StationaryRandomFields
using KernelFunctions: kernelmatrix

grid = FieldGrid((8, 8); spacing = (0.5, 0.5))
power = SinglePowerLaw{2}(-1.2; dc = 0.3)
latent = UnivariateGaussianRandomField(size(grid); σ = 1.4)
covariance = PeriodicCovariance(grid, power, latent)
kernel = spectral_kernel(covariance)
points = [(0.0, 0.0), (0.5, 1.0), (0.13, 0.27)]
K = kernelmatrix(kernel, points)
```

Coordinates can be tuples or vectors. Scalars are also accepted for 1D grids.
For a matrix whose columns contain points, use `kernelmatrix(kernel, X; obsdim=2)`.
The generic `kernelmatrix`, `kernelmatrix!`, and diagonal interfaces are available.
Each pair evaluation sums ``M`` modes, so assembly of an ``m\times m`` matrix costs
``O(m^2M)`` and stores ``O(m^2)`` values. There is no FHT shortcut in this wrapper.
Loading this extension does not automatically accelerate GP inference, linear solves,
log determinants, irregular observations, or observation noise handling.

## Fast action on the full grid

Apply the same covariance to an array on the complete periodic grid with two
orthonormal Hartley transforms:

```julia
x = reshape(sin.(1:64), size(grid))
prepared = prepare(covariance, x)
y = apply(prepared, x)
output = similar(x)
apply!(output, prepared, x)
```

This evaluates ``Q\operatorname{diag}(v\lambda)Q^Tx`` in ``O(M\log M)`` time without
forming the dense matrix. Array order is Julia's column-major `vec` order. Planning
is separate from execution. The current implementation uses temporary arrays; the
buffered API does not promise allocation-free execution.

CPU, CUDA, and Reactant use the same action. Preparation copies frequency metadata
to the input backend once. Updated spectrum parameters can be passed as
`apply(prepared, x, updated_power)`. For pair evaluation, use
`covariance_value(covariance, x, y, updated_power)`. Both paths reevaluate the power
rather than keeping a parameter-dependent cache. Gradients are defined away from
hard spectral cutoff transitions; the existing zero-power convention is unchanged.

## Active latent prior parameters

For Gaussian and Laplace priors, location and scale can be passed directly to
`logpdf`, avoiding construction of a Distributions object with traced scalar fields:

```julia
using Distributions
prior = UnivariateGaussianRandomField((8, 8))
z = zeros(8, 8)
location, logscale = 0.2, -0.3
value = logpdf(prior, z; location, scale = exp(logscale), check_args = false)
```

The scale is standard deviation for Gaussian priors and the Laplace scale ``b`` for
Laplace priors. Normal validation remains enabled by default; compiled objectives
can disable it after enforcing a valid parameterization. This evaluates the latent
density. It does not introduce a signal-space density or Jacobian term.
