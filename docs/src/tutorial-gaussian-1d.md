# Gaussian time series

Think of the first array as independent controls for frequency modes. A power law
sets their relative strengths. Multiplying by **the square root of power** changes
those controls into the desired mode amplitudes. The orthonormal Hartley transform
then turns them into samples along a time or position axis.

```@example gaussian_1d
using StationaryRandomFields, Distributions, Random

rng = MersenneTwister(12)

latent = UnivariateGaussianRandomField(Float64, (256,))
z = rand(rng, latent)

# Use unit sample spacing; frequencies are then measured in cycles per sample.
grid = FieldGrid(size(latent); spacing = (1.0,))

# Choose P(f) = (f / 0.05)^(-2) away from DC; the default DC power is zero.
power = 0.05^2 * SinglePowerLaw{1}(-2.0)

transformer = TransformChain(
    SpectralScaling(power),
    HartleyTransform(),
)
transform = prepare(transformer, latent; grid)
signal = apply(transform, z)

# The prior describes the original latent coefficients, before either operation.
logprior = logpdf(latent, z)

nothing # hide
```

`TransformChain` bundles transformations in left-to-right execution order: spectral
amplitude scaling, then Hartley synthesis. Evaluate `transformer` to inspect its
models and settings. `prepare(transformer, latent; grid)` uses the latent shape and element type, with CPU storage by default.
`logpdf` belongs to `latent`, not to the prepared transform.

`FieldGrid` holds sample counts and coordinate spacing. Here spacing is one, so
frequencies are in cycles per sample. With `spacing=(0.1,)` in seconds, 256 samples
span a periodic interval of 25.6 seconds and frequency bins are separated by
1/25.6 Hz. Omit `grid` for unit spacing, or use
`prepare(transformer, latent; spacing=(0.1,))` as a shortcut.
A `prototype` array can select CUDA or Reactant storage without using its values.

Plot the result separately:

```@example gaussian_1d
using CairoMakie

CairoMakie.activate!(; type = "svg")
set_theme!(Theme(; fontsize = 16))

# Plot positive frequencies only, excluding DC and the even-grid Nyquist bin.
frequency = (1:127) ./ 256
power_values = [power_point(power, f) for f in frequency]

fig = Figure(; size = (1000, 650))

ax = Axis(
    fig[1, 1]; title = "1. Independent Gaussian coefficients",
    xlabel = "Packed frequency-mode index", ylabel = "Coefficient"
)
lines!(ax, eachindex(z), z; color = :steelblue)

ax = Axis(
    fig[1, 2]; title = "2. Chosen power spectrum",
    xlabel = "Frequency [cycles / sample]", ylabel = "Discrete power",
    xscale = log10, yscale = log10
)
lines!(ax, frequency, power_values; color = :darkorange, linewidth = 2)

scaled = apply(prepare(SpectralScaling(power), latent; grid), z)

# This array is still in frequency space, after amplitude scaling.
ax = Axis(
    fig[2, 1]; title = "3. Coefficients after spectral scaling",
    xlabel = "Packed frequency-mode index", ylabel = "Coefficient"
)
lines!(ax, eachindex(scaled), scaled; color = :steelblue)

ax = Axis(
    fig[2, 2]; title = "4. Correlated signal after FHT",
    xlabel = "Time [sample units]", ylabel = "Signal [arbitrary units]"
)
lines!(ax, 0:255, signal; color = :purple, linewidth = 2)
save("gaussian-1d.svg", fig)
nothing # hide
```

![Four stages of one-dimensional Gaussian spectral synthesis](gaussian-1d.svg)

The first and third panels show **frequency coefficients**, not light curves. The
last panel is the signal. Power follows an inverse-square law in frequency here, so slow changes
are stronger than rapid ones. Zero DC power removes the constant mode. The plot
shows the target power, not the periodogram of this single draw; individual draws
fluctuate around the target in power. The reference frequency `0.05` sets an amplitude
convention for this example, rather than a universal physical normalization.

The finite grid is periodic: the signal wraps around at its boundary. When using
these samples as a light-curve noise model, choose the grid extent, cadence, and
spectral cutoffs for your observation. See [Timmer & König (1995)](https://scixplorer.org/abs/1995A%26A...300..707T/abstract)
for the original time-series motivation and [Mathematical Basis](math.md) for the
package's discrete-power convention.
