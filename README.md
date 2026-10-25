# StationaryRandomFields.jl
| Documentation                                                                                                                                                                                                                                   | Code Status                                                                                                                                                                                                                                                                                                                                                                              | Guidelines                                                                                                                                                                                                                                                                                                                                               |
| :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://EHTJulia.github.io/StationaryRandomFields.jl/stable/)[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://EHTJulia.github.io/StationaryRandomFields.jl/dev/) | [![Build Status](https://github.com/EHTJulia/StationaryRandomFields.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/EHTJulia/StationaryRandomFields.jl/actions/workflows/CI.yml?query=branch%3Amain)[![Coverage](https://codecov.io/gh/EHTJulia/StationaryRandomFields.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/EHTJulia/StationaryRandomFields.jl) | [![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac) |

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

## Generation of correlated fields

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
