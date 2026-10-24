# Tutorial

## Gaussian correlated signal

```julia
using StationaryRandomFields, Distributions, Random

latent = UnivariateGaussianRandomField(Float64, (32, 32))
z = rand(MersenneTwister(12), latent)
grid = FieldGrid(size(z); spacing = (0.1, 0.1))
ps = SaturatedPowerLaw{2}(-2.0, 0.1, 4.0)
spectral = prepare(SpectralScaling(grid, ps), z)
hartley = prepare(HartleyTransform(size(z)), z)
signal = apply(hartley, apply(spectral, z))
logprior = logpdf(latent, z)
```

Update a spectrum without rebuilding a plan or caching stale amplitudes:

```julia
updated = 2.0 * ps + ConstantPowerSpectrum{2}(0.01)
signal = apply(hartley, apply(spectral, z, updated))
```

## Laplace coefficients

```julia
latent = UnivariateLaplaceRandomField(Float32, (32, 32); b = 0.2f0)
z = rand(latent)
synthesis = prepare(WaveletTransform(size(z); levels = 3), z)
signal = apply(synthesis, z)
logprior = logpdf(latent, z)
```

Use `IdentityTransform(size(z))` for a direct image-space Laplace prior. Use
`TransformChain(ScaleTransform(size(z), weights), WaveletTransform(size(z); levels=3))`
for coefficient weights. The weights must broadcast to the coefficient shape.

## Sampling a transformed field

```julia
field = TransformedRandomField(latent, synthesis)
x = rand(field)
```

`logpdf(field, x)` is deliberately unavailable. Evaluate `logpdf(field.latent, z)`.
For a vector-valued field, `rand(field, n)` follows Distributions' matrix-of-columns
convention; higher-rank fields use arrays of sample arrays. `rand!(rng, field, batch)`
also accepts a numeric array whose trailing axes index independent samples.

## Differentiation and compiled execution

Prepare transforms outside the objective and keep random draws fixed during AD.
`examples/latent_inference.jl` contains a complete latent log-posterior example.
The CPU Enzyme suite checks latent, wavelet, Hartley, and spectral-parameter gradients.

For Reactant, transfer the latent array, prepare using that array as the prototype, and
compile `apply` or the objective. Trainable parameters should be device arrays, not Julia
scalar constants. A scalar parameter can be represented by a one-element array and
reduced with `sum`, avoiding scalar indexing. When constructing a dynamic positive
constant spectrum with `exp(sum(parameters))`, use `check_args=false` to keep host
validation outside the compiled objective. Validate physical parameter ranges before
compilation and use positive parameterizations for optimization.

CPU and CUDA sampling use backend-specific RNGs. Compiled Reactant sampling uses an
explicit `ReactantRNG`, preserving state between calls. RNG streams and sample values
need not match between backends; numerical transforms are compared with fixed inputs.

For trainable spectrum geometry, construct modifiers with
`Stretch(exp(logwidth1), exp(logwidth2); check_args=false)` and `Rotate(angle)`,
then combine them with `modify`. The compiled gradient suite checks both widths and
the two-dimensional rotation angle against finite differences. Static proper rotation
matrices support arbitrary spectrum rank.
