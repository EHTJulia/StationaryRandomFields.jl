# Optional visualization of the fitted FITS image produced by closure_imaging.jl.
using Comrade, CairoMakie

length(ARGS) in (1, 2) || error("Usage: julia --project=. plot_reconstruction.jl input.fits [output.png]")

# Read the saved image and place the plot beside it unless a path is supplied.
input = ARGS[1]
output = length(ARGS) == 2 ? ARGS[2] : joinpath(dirname(input), "wavelet_closure_map.png")
image = load_fits(input, IntensityMap)

# The optional tutorial style uses Comrade's image visualization defaults.
if get(ENV, "WAVELET_TUTORIAL_STYLE", "0") == "1"
    figure = imageviz(image; size = (600, 500))
    save(output, figure)
    println("Comrade tutorial-style plot saved: ", output)
else
    # Display the same FITS pixels with explicit axes and a flux color bar.
    figure = Figure(; size = (700, 600))
    axis = Axis(
        figure[1, 1]; xlabel = "Relative right ascension (μas)", ylabel = "Relative declination (μas)",
        title = get(ENV, "WAVELET_PLOT_TITLE", "Laplace-wavelet closure MAP"), xreversed = true, aspect = DataAspect(),
    )
    map = heatmap!(axis, rad2μas.(collect(image.X)), rad2μas.(collect(image.Y)), baseimage(image); colormap = :afmhot)
    Colorbar(figure[1, 2], map; label = "Jy / pixel")
    save(output, figure)
    println("Plot saved: ", output)
end
