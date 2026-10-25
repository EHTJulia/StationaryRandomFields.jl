# Gaussian image textures

Compare a single power law, a saturated power law, Matérn power, and Markov power
using the same Gaussian coefficients. Each spectrum is normalized to unit expected
pointwise variance on the same periodic grid.

```@example gaussian_images
using StationaryRandomFields, Random

latent = UnivariateGaussianRandomField((128, 128))
z = rand(MersenneTwister(23), latent)

models = (
    SinglePowerLaw{2}(-3.0),
    SaturatedPowerLaw{2}(-3.0, 0.035, 0.18),
    MaternPS{2}(8.0, 1.5),
    MarkovPS{2}((8.0,)),
)

images = map(models) do model
    transformer = TransformChain(
        SpectralScaling(model; normalization = :unit_variance),
        HartleyTransform(),
    )
    transform = prepare(transformer, latent)
    apply(transform, z)
end
nothing # hide
```

`TransformChain` applies its transformations from left to right: spectral scaling,
then Hartley synthesis. Each model supplies its own prepared chain.

The shared latent draw and color scale make the texture differences visible.
Normalization matches ensemble variance, not the sample variance of each image.

```@example gaussian_images
using CairoMakie

CairoMakie.activate!()
set_theme!(Theme(; fontsize = 16))
titles = ("Single power law", "Saturated power law", "Matérn", "Markov")
limit = maximum(img -> maximum(abs, img), images)
fig = Figure(; size = (900, 800))
for (j, (img, title)) in enumerate(zip(images, titles))
    local row, col = fld(j - 1, 2) + 1, mod1(j, 2)
    local ax = Axis(fig[row, col]; title, aspect = DataAspect())
    heatmap!(ax, img; colormap = :balance, colorrange = (-limit, limit))
end
Colorbar(fig[:, 3]; colormap = :balance, limits = (-limit, limit), label = "Signal")
save("gaussian-textures.png", fig)
nothing # hide
```

![Four Gaussian image textures with shared latent coefficients and color scale](gaussian-textures.png)

- **Single power law:** concentrates power in broad structures.
- **Saturated power law:** flattens below the inner frequency and removes power above the outer frequency.
- **Matérn:** controls correlation range and smoothness separately.
- **Markov:** the one-term spectrum retains more fine-scale texture.

These examples compare particular parameter choices. See [Spectral transforms](math-spectral.md)
for the formulas, frequency units, and DC conventions. Change the models and rerun the synthesis block to compare other textures.
