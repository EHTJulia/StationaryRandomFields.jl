# Daubechies wavelet signals

Now coefficients represent features at different locations and scales. A Laplace
prior prefers small coefficients, with heavier tails than a Gaussian of the same
variance. Wavelet synthesis assembles these features into the signal. Here we use
the periodic Daubechies 4 (Db4) wavelet with the maximum decomposition depth.

```@example wavelet_1d
using StationaryRandomFields, Distributions, Random, Wavelets

# b is the Laplace scale, not its standard deviation.
latent = UnivariateLaplaceRandomField(Float64, (256,); b = 0.2)
z = rand(MersenneTwister(34), latent)

# Use all eight dyadic levels for 256 samples.
synthesis = prepare(WaveletTransform(wavelet(WT.db4)), latent)

signal = apply(synthesis, z)

logprior = logpdf(latent, z)

nothing # hide
```

Plot the result separately:

```@example wavelet_1d
using CairoMakie

CairoMakie.activate!(; type = "svg")
set_theme!(Theme(; fontsize = 16))

fig = Figure(; size = (1000, 370))
ax = Axis(
    fig[1, 1]; title = "Independent Laplace coefficients",
    xlabel = "Packed wavelet-coefficient index", ylabel = "Coefficient"
)
stem!(ax, eachindex(z), z; color = :steelblue, markersize = 3)

ax = Axis(
    fig[1, 2]; title = "Signal after Db4 wavelet synthesis",
    xlabel = "Position [sample units]", ylabel = "Signal [arbitrary units]"
)
lines!(ax, 0:255, signal; color = :purple)
save("laplace-1d.svg", fig)
nothing # hide
```

![Laplace coefficients and their synthesized one-dimensional signal](laplace-1d.svg)

Coefficients are packed by wavelet scale. See [Wavelet transforms](math-wavelet.md)
for the packing convention and the role of the Laplace prior.
