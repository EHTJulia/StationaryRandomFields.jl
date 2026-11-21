# KernelFunctions.jl extension

Load `KernelFunctions` to use the covariance of a Gaussian spectral field as a kernel.
The kernel is periodic; its definition is in [Spectral transforms](math-spectral.md#Periodic-covariance).

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
``O(m^2M)`` and stores ``O(m^2)`` values.

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
forming the dense matrix. Array order follows Julia's `vec` convention.

Pass `updated_power` to `apply(prepared, x, updated_power)` to change the spectrum
without rebuilding the transform plans.
