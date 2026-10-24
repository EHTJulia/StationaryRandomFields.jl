# Introduction

Astronomical signals often have structure across neighboring samples. A variable source
can remain bright for several time samples; diffuse emission changes gradually across
an image; a turbulent screen has fluctuations on many spatial scales. Drawing a new,
independent noise value at every pixel misses that structure.

StationaryRandomFields.jl generates **correlated noise for signals of any array dimension
with specified power spectra**, and **signals synthesized from a sparsity-promoting
wavelet space**. It also provides these generative models as building blocks for Bayesian
reconstruction. The same machinery can describe noise added to a signal or a prior on
the unknown signal itself.

## Two useful ways to build a signal

| What you want to describe | Start with | Then apply |
|---|---|---|
| Correlated fluctuations with a chosen distribution of power over scales | Independent Gaussian coefficients | Multiply by the square root of the power spectrum, then a Hartley transform |
| A signal described economically by wavelet coefficients | Independent Laplace coefficients | Optional coefficient weights, then wavelet synthesis |

A **latent field** is simply the array of random numbers before these operations. A
**signal field** is the result you would plot, pass to a telescope model, or compare with
data. In the Gaussian route, the latent array indexes frequency modes. In the wavelet
route, it stores coefficients at different locations and scales. Neither array should
be mistaken for an observed image before synthesis.

A power spectrum describes how much variation to expect at each frequency: low
frequencies make broad, slowly varying structures, while high frequencies make fine
texture or rapid changes. A negative power-law index favors large scales. Wavelets
instead describe localized features at several scales. A Laplace prior favors many
small coefficients while allowing occasional larger ones; a random draw is not exactly
zero at most entries, and is not a thresholded sparse reconstruction.

## Why this is useful in astronomy

- **Time-domain simulations:** generate correlated variability or noise realizations to
  test a light-curve analysis. The spectral-synthesis motivation goes back to
  [Timmer & König (1995)](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..707T/abstract).
- **Radio imaging priors:** describe extended emission or fluctuations around a reference
  image, and infer the correlation parameters along with the image. Related Bayesian
  imaging work includes [Arras et al. (2021)](https://doi.org/10.1051/0004-6361/202039258)
  and [Tiede et al. (2026)](https://arxiv.org/abs/2511.17706). The latter uses Markov-field
  image priors; it illustrates the imaging motivation rather than an identical spectral model.
- **Scattering and atmospheric phase screens:** use a chosen spectrum to draw a screen
  of phase fluctuations, for example as an ingredient in Sgr A* scattering simulations
  ([Johnson 2016](https://arxiv.org/abs/1610.05326)) or empirical interferometric
  atmospheric/ionospheric models ([Chege et al. 2021](https://doi.org/10.1017/pasa.2021.22)).
  Propagation, wavelength scaling, calibration, and the visibility response belong in
  the surrounding physical model.
- **Sparse Bayesian imaging:** place a Laplace prior on wavelet coefficients to favor
  images with a compact multiscale description. The corresponding negative log prior
  is an L1 penalty in coefficient space. This supplies a prior, rather than a complete
  reconstruction algorithm.

The default Gaussian signal is signed. To model positive sky brightness, choose an
additional physical mapping, such as exponentiation, in your observation model.
The plotting examples use arbitrary units; select spacings, cutoffs, amplitudes, and
units appropriate to the science before interpreting a field as a physical screen.

## From a model to inference

`UnivariateRandomField` stores an array shape and the shared scalar distribution of
its independent latent coefficients. An `AbstractFieldTransform` describes how those
coefficients become a signal. `TransformedRandomField` combines them for sampling.
The output shape may differ from the latent shape for a custom transform.

In Bayesian use, estimate the latent coefficients and transform parameters, evaluate
`logpdf` on the latent field, and evaluate the likelihood on the synthesized signal.
Signal-space `logpdf` is not exposed. See the [tutorial](tutorial.md) for working examples
and the [mathematical basis](math.md) for the explicit probability model.

For Gaussian spectral models, the optional [KernelFunctions interface](kernels.md)
provides the covariance of the same finite periodic-grid sampler.

## Preparing work for a backend

`prepare(transform, prototype)` builds reusable plans and static metadata for an array
with the prototype's type and shape. Prepare outside differentiated or compiled
objectives; model-dependent spectral amplitudes are recomputed during execution.
Array storage selects CPU, CUDA, or Reactant execution. Transfer coefficient weights
explicitly to the input's backend. Prepare separate objects and output buffers for
concurrent tasks; Hartley `apply!` rejects overlapping input and output storage.

`rand(field)` follows the Distributions CPU-allocation convention. For device-prepared
fields, use `rand!(rng, field, device_output)`. CPU sampling accepts an `AbstractRNG`,
CUDA uses its own RNG, and compiled Reactant sampling uses an explicit `ReactantRNG`.
Keep RNG state outside differentiation. The [tutorial](tutorial.md) explains how to
hold a random draw fixed while differentiating a model.
