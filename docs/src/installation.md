# Installation

Use Julia 1.10 LTS or newer. Download Julia from [julialang.org](https://julialang.org/downloads/).

## Registered release

In a Julia environment for your project, install the registered release with:

```julia
using Pkg
Pkg.add("StationaryRandomFields")
```

These pages describe the redesigned development API. The registered 0.1 series has the
older noise-generator API; its documentation is available at
[version 0.1.3](https://ehtjulia.github.io/StationaryRandomFields.jl/v0.1.3/).
Installing that release does not provide the new examples on this development site.

For the tutorial plots, install `CairoMakie` separately. CUDA, Reactant, and
KernelFunctions are optional integrations; install and load the library you intend to use.
