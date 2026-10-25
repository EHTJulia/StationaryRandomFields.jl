# Daubechies wavelet images

For an image, we can additionally give broad features more amplitude than fine
texture. `ScaleTransform` weights coefficients before wavelet synthesis. In the
periodic packing used here, the approximation coefficients occupy the low-index
corner at each level; the loop assigns larger weights to successive coarse corners.
We use periodic Db4 and all seven available levels. An explicit `levels` keyword
can limit the depth; the default is the deepest common dyadic level supported by all axes.

```@example wavelet_2d
using StationaryRandomFields, Distributions, Random, Wavelets

latent = UnivariateLaplaceRandomField(Float64, (128, 128); b = 0.2)
z = rand(MersenneTwister(45), latent)

# Begin with small amplitudes at all scales, then strengthen coarse-scale blocks.
dwt = WaveletTransform(wavelet(WT.db4))
levels = prepare(dwt, latent).levels # Seven levels for a 128-by-128 image.
weights = fill(0.1, size(z))
# Each successive low-index corner contains the next coarser approximation.
for level in 1:levels
    width = size(z, 1) ÷ 2^level
    weights[1:width, 1:width] .= 0.1 * 2.0^level
end

transform = TransformChain(
    ScaleTransform(weights), dwt
)
synthesis = prepare(transform, latent)
signal = apply(synthesis, z)

# Weighting is part of the transform; the prior still belongs to the unweighted z.
logprior = logpdf(latent, z)

nothing # hide
```

Plot the result separately:

```@example wavelet_2d
using CairoMakie

CairoMakie.activate!(; type = "svg")
set_theme!(Theme(; fontsize = 16))

function field_panel!(slot, values, title, xlabel, ylabel; colormap = :balance)
    ax = Axis(slot[1, 1]; title, xlabel, ylabel, aspect = DataAspect())
    image = heatmap!(ax, values; colormap)
    Colorbar(slot[1, 2], image)
    return ax
end

fig = Figure(; size = (1200, 390))
field_panel!(fig[1, 1], z, "1. Laplace latent coefficients", "Packed index 1", "Packed index 2")
field_panel!(fig[1, 2], weights .* z, "2. Coefficients with scale weights", "Packed index 1", "Packed index 2")
field_panel!(fig[1, 3], signal, "3. Db4-synthesized signal", "Pixel x", "Pixel y")
save("laplace-2d.png", fig)
nothing # hide
```

![Laplace coefficients, scale weighting, and a two-dimensional wavelet signal](laplace-2d.png)

Coarse-scale weights favor broad structure. During inference, evaluate the prior
on the unweighted latent coefficients; see [Latent fields](math-latent.md).
