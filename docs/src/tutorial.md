# Tutorials

Each example is self-contained. Install `CairoMakie` for the plots.

- [Gaussian time series](tutorial-gaussian-1d.md)
- [Gaussian image textures](tutorial-gaussian-2d.md): compare four spectral shapes.
- [Daubechies wavelet signals](tutorial-wavelet-1d.md)
- [Daubechies wavelet images](tutorial-wavelet-2d.md)
- [Custom fields](tutorial-custom-field.md)
- [Custom transforms](tutorial-custom-transform.md)

For sampling and prior evaluation, see [Latent fields and inference](math-latent.md).

Transform descriptors specify operations without a data shape. Prepare them with
`prepare(transformer, latent)`, then use `apply(transform, z)` on a latent draw.
The [first Gaussian tutorial](tutorial-gaussian-1d.md) explains `FieldGrid`, spacing,
and backend prototypes. Only `logpdf(latent, z)` evaluates the prior; transforms
do not define a probability density.
