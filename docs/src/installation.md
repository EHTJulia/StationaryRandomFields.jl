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

## This development checkout

Until the redesigned package and companion Hartley packages are released, use the
reviewed development checkout. From its repository root, run:

```sh
julia tools/setup.jl
julia --project -e 'using StationaryRandomFields'
```

The setup script installs the companion commits recorded in `tools/setup.jl` and
instantiates the package environment. It does not install plotting dependencies into
the runtime package. If you already maintain sibling Hartley checkouts, see the
[repository README](https://github.com/EHTJulia/StationaryRandomFields.jl#development-and-validation)
for the local `Pkg.develop` alternative.

## Run the figures and build this site

The documentation environment includes CairoMakie for headless plotting. From the
same repository root:

```sh
julia tools/setup.jl docs
julia --project=docs docs/make.jl
```

This executes the tutorial examples and builds the VitePress site in
`docs/build/final_site`. To preview it locally with its configured site prefix:

```sh
cd docs
npm run docs:preview
```

Open the URL printed by VitePress. The site uses `/StationaryRandomFields.jl/` as its
local base path. No deployment occurs unless `SRF_DEPLOY_DOCS=true` is set explicitly.
For plotting in your own project, add `CairoMakie` separately with `Pkg.add("CairoMakie")`.

CUDA, Reactant, and KernelFunctions are optional integrations: install and load the
backend or library you intend to use. The basic CPU examples require none of them.
