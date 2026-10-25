# Latent fields and inference

A latent field is an array of random coefficients ``z``. A transform maps it to
a signal ``x=T_\theta(z)``. The prior density is evaluated on ``z``.

## Independent coefficients

For `M = prod(size(latent))` independent elements:

```math
\log p_Z(z)=\sum_{i=1}^M\log p(z_i).
```

For a Laplace marginal with location ``\mu`` and scale ``b``:

```math
\log p_Z(z)=-M\log(2b)-\frac{1}{b}\sum_i|z_i-\mu|.
```

The scale is not a standard deviation: its variance is ``2b^2``. At the Laplace cusp,
AD uses the zero subgradient.

## Inference

```math
\log p(z,\theta\mid y)=\log p(y\mid T_\theta(z))+
\log p_Z(z)+\log p(\theta)+C.
```

Infer the latent coefficients ``z`` and model parameters ``\theta``. Evaluate
`logpdf(latent, z)` in latent space and the likelihood on the transformed signal.

## Sampling

```@example inference
using StationaryRandomFields, Distributions, Wavelets

latent = UnivariateLaplaceRandomField((128, 128); b = 0.2)
prototype = zeros(size(latent))
synthesis = prepare(WaveletTransform(wavelet(WT.db4)), latent; prototype)

field = TransformedRandomField(latent, synthesis)

x = rand(field)
size(x)
```

Use `logpdf(field.latent, z)` to evaluate the prior of the coefficients.

## Location and scale

Pass location and scale directly when fitting prior parameters:

```@example inference
using Distributions
prior = UnivariateGaussianRandomField((8, 8))
z = zeros(8, 8)
location, logscale = 0.2, -0.3
value = logpdf(prior, z; location, scale = exp(logscale))
```

`scale` is the Gaussian standard deviation or the Laplace scale ``b``.
Use a positive parameterization such as `exp(logscale)` when fitting it.

Prepare transforms outside the objective and keep latent draws fixed while differentiating.
See [custom fields](tutorial-custom-field.md) for other latent distributions.
