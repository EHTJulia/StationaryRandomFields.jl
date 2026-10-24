```@raw html
---
layout: home
hero:
  name: "StationaryRandomFields.jl"
  text: "From random coefficients to structured signals"
  tagline: "Correlated noise with specified power spectra, and signals synthesized from sparse wavelet priors — in any array dimension."
  actions:
    - theme: brand
      text: Get started
      link: /installation
    - theme: alt
      text: See the illustrated tutorial
      link: /tutorial
    - theme: alt
      text: GitHub
      link: https://github.com/EHTJulia/StationaryRandomFields.jl
features:
  - title: Correlated signals
    details: Turn independent Gaussian coefficients into time series, images, and higher-dimensional fields with a chosen power spectrum.
  - title: Sparse wavelet priors
    details: Synthesize signals from Laplace-distributed coefficients, with optional weights for different scales.
  - title: Bayesian building blocks
    details: Evaluate priors in latent space, differentiate model parameters, and connect finite-grid covariances to KernelFunctions.jl.
---
```

## What can I make with it?

Generate realistic correlated noise for signal data following specified power spectra,
or build signals from a wavelet space with a sparsity-promoting prior. A one-dimensional
array can represent a light curve; a two-dimensional array can represent an image or a
phase screen; higher-dimensional arrays can represent data cubes.

Start with the [motivation and astronomical applications](introduction.md),
[install the package](installation.md), and follow the [illustrated tutorial](tutorial.md).
The [mathematical basis](math.md) explains the normalization and probability model;
[covariance kernels](kernels.md) and the [API reference](autodocs.md) cover further use.
