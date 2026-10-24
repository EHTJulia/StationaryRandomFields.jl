# StationaryRandomFields.jl

Generate realistic correlated noise for signal data following specified power spectra,
and synthesize signals from a sparsity-promoting wavelet space, in any array dimension.
Use a one-dimensional array for a light curve, two dimensions for an image or phase
screen, or higher dimensions for data cubes. These same models provide latent priors
for Bayesian signal reconstruction.

Gaussian coefficients are scaled by the square root of a power spectrum and mapped
into signal space with a multidimensional Hartley transform. Laplace coefficients
provide a sparsity-promoting wavelet prior, optionally weighted by scale. Applications
include correlated variability, radio imaging priors, scattering/atmospheric phase
screens, and sparse Bayesian imaging. Physical propagation and telescope responses
are supplied by the application.

## Installation

For the registered release:

```julia
using Pkg
Pkg.add("StationaryRandomFields")
```

The examples below describe the redesigned development API, which differs from the
registered 0.1 series. For this development checkout, run `julia tools/setup.jl` from
the repository root, then start Julia with `julia --project`. The script installs the
required companion revisions. See [installation](docs/src/installation.md) for details
and [the illustrated tutorial](docs/src/tutorial.md) for executable CairoMakie figures.

## A Gaussian correlated field

```julia
using StationaryRandomFields, Distributions, Random

latent = UnivariateGaussianRandomField(Float32, (64, 64))
z = rand(MersenneTwister(12), latent)
ps = SaturatedPowerLaw{2}(-2.0f0, 0.1f0, 4.0f0)
transform = TransformChain(SpectralScaling(ps), HartleyTransform())
prepared = prepare(transform, latent; spacing = (0.1, 0.1))
x = apply(prepared, z)
logprior = logpdf(latent, z)
field = TransformedRandomField(latent, prepared)
sample = rand(field)
```

`UnivariateRandomField` describes the common scalar distribution of independent latent
elements. A sample is an N-dimensional array; the full prior is their joint product.
`logpdf` is evaluated **only in latent space**. Signal-space `logpdf` and general signal
support membership deliberately throw informative errors. No Jacobian API is required.

For a wavelet-coefficient Laplace prior with all available dyadic levels:

```julia
latent = UnivariateLaplaceRandomField(Float32, (64, 64); b = 0.2f0)
z = rand(latent)
synthesis = prepare(WaveletTransform(), latent)
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
Pkg.develop(
    [
        PackageSpec(path = "../AbstractFastHartleyTransforms"),
        PackageSpec(path = "../FastHartleyTransform"),
    ]
)
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

[Runic](https://github.com/fredrikekre/Runic.jl) is pinned to 1.13.0 in the tooling
environment. The check is read-only and covers Julia files, docstrings, and Julia code
blocks in Markdown, including Documenter examples. Code comments,
docstrings, and configuration comments use American English. Contributions follow
[ColPrac](https://github.com/SciML/ColPrac), with topic branches and independent review.


The optional KernelFunctions extension wraps `PeriodicCovariance` as a finite periodic
spectral kernel. Dense kernel matrices and the fast full-grid FHT covariance action
are separate APIs. See [covariance kernels](docs/src/kernels.md) for normalization,
off-grid semantics, parameter differentiation, and examples.
