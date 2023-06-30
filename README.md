# StationaryRandomFields.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://EHTJulia.github.io/StationaryRandomFields.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://EHTJulia.github.io/StationaryRandomFields.jl/dev/)
[![Build Status](https://github.com/EHTJulia/StationaryRandomFields.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/EHTJulia/StationaryRandomFields.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/EHTJulia/StationaryRandomFields.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/EHTJulia/StationaryRandomFields.jl)


This Julia package is designed to provide toolkits for handling signals sampled from stationary Gaussian random fields that follow a specific stationary power spectrum in the Fourier domain. The package offers:

- Abstract types and methods for multi-dimensional discrete signals.
- Abstract multi-dimensional power spectrum types, which can be easily scaled, stretched, and rotated.
- toolkits to sample signals from input (or generated) Gaussian random fields in the Fourier space.

The package is currently in the significant stages of development. We are planning to implement:
- Chain rules to make the relevant methods compatible with Julia's standard automatic differentiation packages.
- Custom kernels for [KernelFunctions.jl](https://github.com/JuliaGaussianProcesses/KernelFunctions.jl), enabling the integration of stationary Gaussian random fields characterized by power spectra into the rich ecosystem of the [Julia Gaussian Processing packages](https://github.com/JuliaGaussianProcesses).

**Disclaimer**: The package may have significant changes in near future given the early stage of the development. 


## Installation
The package is being registered to the Julia standard repository as of June 30, 2023. Once registered, it will be installable by
```julia
using Pkg
Pkg.add("StationaryRandomFields")
```

## Usage
To be written soon.


## Acknowledgements
The development & developers of this package has been finantially supported by the following programs.
- AST-2107681, National Science Foundation, USA: v0.1.0 - present
- AST-2034306, National Science Foundation, USA: v0.1.0 - present
- OMA-2029670, National Science Foundation, USA: v0.1.0 - present
- AST-1950348, National Science Foundation, USA: v0.1.0 - present
- AST-1935980, National Science Foundation, USA: v0.1.0 - present