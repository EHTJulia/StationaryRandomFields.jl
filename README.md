# StationaryRandomFields.jl

Array-valued latent priors and composable transforms for Bayesian signal reconstruction.
The current redesign provides independent Gaussian/Laplace latent fields, power spectra,
orthonormal multidimensional Hartley synthesis, and multilevel wavelet synthesis.

```julia
using StationaryRandomFields, Distributions, Random

latent = UnivariateGaussianRandomField(Float32, (64, 64))
z = rand(MersenneTwister(12), latent)
grid = FieldGrid(size(z); spacing=(0.1, 0.1))
ps = SaturatedPowerLaw{2}(-2f0, 0.1f0, 4f0)
transform = TransformChain(SpectralScaling(grid, ps), HartleyTransform(size(z)))
prepared = prepare(transform, z)
x = apply(prepared, z)
logprior = logpdf(latent, z)
field = TransformedRandomField(latent, prepared)
sample = rand(field)
```

`UnivariateRandomField` describes the common scalar distribution of independent latent
elements. A sample is an N-dimensional array; the full prior is their joint product.
`logpdf` is evaluated **only in latent space**. Signal-space `logpdf` and general signal
support membership deliberately throw informative errors. No Jacobian API is required.

For a wavelet-coefficient Laplace prior:

```julia
latent = UnivariateLaplaceRandomField(Float32, (64, 64); b=0.2f0)
z = rand(latent)
synthesis = prepare(WaveletTransform(size(z); levels=3), z)
x = apply(synthesis, z)
logprior = logpdf(latent, z)
```

`TransformChain(a, b)` applies `a` first and `b` second. `IdentityTransform` expresses
an image-space prior; `ScaleTransform` expresses scalar or coefficient-wise amplitude
weights. These weights are distinct from a Fourier power spectrum.

See [the tutorial](docs/src/tutorial.md), [the mathematical conventions](docs/src/math.md),
[the redesign review](REVIEW.md), and [the reproducible measurements](benchmark/RESULTS.md).
An executable inference example is in `examples/latent_inference.jl`.

## Development and validation

Julia 1.10 LTS or newer is required. The companion Hartley packages currently need the
reviewed development commits listed in `REVIEW.md`; they are not replaced by an old
registry release. The checked-in development manifest uses sibling checkouts. With the
three repositories checked out next to each other, run:

```julia
using Pkg
Pkg.develop([
    PackageSpec(path="../AbstractFastHartleyTransforms"),
    PackageSpec(path="../FastHartleyTransform"),
])
Pkg.instantiate()
Pkg.test()
```

For standalone CI, `tools/setup.jl` uses explicit companion commit revisions.
CPU tests cover Float32/Float64, numerical references, shape validation, and Enzyme.
The Reactant suite runs compiled execution and differentiation on CPU by default.
CUDA and Reactant GPU suites require an actual device and are run explicitly; an
unavailable device is an error, not a silently passing test. No DMOG job is launched
by the repository scripts or workflows.

```sh
julia --project=tools -e 'using Pkg; Pkg.instantiate()'
julia --project=tools tools/format.jl --write
julia --project=tools tools/format.jl
```

JuliaFormatter is pinned to 2.1.2 and uses the root BlueStyle configuration. Code comments,
docstrings, and configuration comments use American English. Contributions follow
[ColPrac](https://github.com/SciML/ColPrac), with topic branches and independent review.
