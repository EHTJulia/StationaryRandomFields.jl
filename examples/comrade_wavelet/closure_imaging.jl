# Closure-only imaging with Laplace coefficients and a Sym6 wavelet prior.
using Comrade, VLBIFiles, VLBIImagePriors, Distributions, StationaryRandomFields
using Wavelets, Enzyme, Optim, LogDensityProblems, Random

# `l` is the coefficient array itself, not the `c.params` field used by the
# GMRF prior in Comrade's original closure-imaging tutorial.
@sky function wavelet_sky(grid; mean_image, coefficient_prior, synthesis)
    l ~ coefficient_prior
    σimg ~ VLBIExponential(0.1)
    fg ~ VLBIUniform(0.0, 1.0)

    # Synthesize a fluctuation and apply it to the mean image as a log ratio.
    coefficients = reshape(l, size(mean_image))
    fluctuation = apply(synthesis, coefficients)
    raster = apply_fluctuations(CenteredLR(), mean_image, σimg .* fluctuation)

    # The compact and extended components have fractions 1-fg and fg.
    compact = ContinuousImage((1 - fg) .* raster, BSplinePulse{3}())
    extended = modify(Gaussian(), Stretch(μas2rad(250.0), μas2rad(250.0)), Renormalize(fg))
    return compact + extended
end

function main(; mean_image = nothing)
    started = time()

    # Read the optimization settings and output location.
    seed = parse(Int, get(ENV, "WAVELET_SEED", "20261031"))
    rng = MersenneTwister(seed)
    npix = parse(Int, get(ENV, "WAVELET_NPIX", "32"))
    maxiters = parse(Int, get(ENV, "WAVELET_MAXITERS", "30"))
    ispow2(npix) || throw(ArgumentError("WAVELET_NPIX must be a power of two."))
    output = get(ENV, "WAVELET_OUTPUT", joinpath(@__DIR__, "output"))
    mkpath(output)

    # The Comrade installation may include its public tutorial data. A user can
    # override the path without copying or changing the original data file.
    default_data = joinpath(pkgdir(Comrade), "examples", "Data", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits")
    data_path = get(ENV, "WAVELET_UVFITS", default_data)
    isfile(data_path) || error("UVFITS data not found: $data_path; set WAVELET_UVFITS to the tutorial data path.")
    uvd = VLBIFiles.load(VLBIFiles.UVData, data_path)

    # Form scan-averaged closure quantities and add 2% fractional noise.
    dlcamp, dcphase = extract_table(
        uvd,
        LogClosureAmplitudes(; time_average = VLBI.GapBasedScans()),
        ClosurePhases(; time_average = VLBI.GapBasedScans()),
    )
    add_fractional_noise!(dlcamp, 0.02)
    add_fractional_noise!(dcphase, 0.02)

    # Use a normalized Gaussian mean unless a matching image was supplied.
    grid = imagepixels(μas2rad(150.0), μas2rad(150.0), npix, npix)
    fwhmfac = 2sqrt(2log(2))
    mean_model = modify(Gaussian(), Stretch(μas2rad(50.0) / fwhmfac))
    default_mean = intensitymap(mean_model, grid)
    if isnothing(mean_image)
        mean_image = default_mean
    else
        size(mean_image) == size(default_mean) || throw(ArgumentError("The mean image must match the reconstruction grid size."))
        collect(mean_image.X) ≈ collect(default_mean.X) || throw(ArgumentError("The mean image X coordinates must match the reconstruction grid."))
        collect(mean_image.Y) ≈ collect(default_mean.Y) || throw(ArgumentError("The mean image Y coordinates must match the reconstruction grid."))
        all(isfinite, baseimage(mean_image)) && all(>=(0), baseimage(mean_image)) || throw(ArgumentError("The mean image must contain finite, nonnegative values."))
        Comrade.flux(mean_image) > 0 || throw(ArgumentError("The mean image must have positive flux."))
    end
    mean_image = mean_image ./ Comrade.flux(mean_image)

    # The latent distribution is independent Laplace(0,1), whose variance is 2.
    # `Product` gives Comrade an ordinary array-valued prior/transport; its
    # density is the same as this StationaryRandomFields latent field.
    dims = (npix, npix)
    latent = UnivariateLaplaceRandomField(Float64, dims; b = 1.0)
    coefficient_prior = product_distribution(fill(Laplace(0.0, 1.0), prod(dims)))
    # Preserve the seed sequence used by the documented reconstruction.
    rand(rng, latent)

    # Periodic orthogonal wavelet synthesis; larger low-index coarse blocks are
    # weighted more strongly. The prior remains on the unweighted coefficients.
    dwt = WaveletTransform(dims, wavelet(WT.sym6))
    levels = prepare(dwt, zeros(dims)).levels
    weights = fill(0.25, dims)
    for level in 1:levels
        width = npix ÷ 2^level
        weights[1:width, 1:width] .= 0.25 * 2.0^(level / 2)
    end
    synthesis = prepare(TransformChain(ScaleTransform(dims, weights), dwt), zeros(dims))

    # Build the posterior and map its parameters to flat optimizer coordinates.
    skym = wavelet_sky(grid; mean_image, coefficient_prior, synthesis)
    post = VLBIPosterior(skym, dlcamp, dcphase; admode = Enzyme.set_runtime_activity(Enzyme.Reverse))
    tpost = asflat(post)
    initial = prior_sample(rng, post)
    x0 = Comrade.inverse(tpost, initial)
    objective(x) = -logdensityof(tpost, x)
    function gradient!(g, x)
        _, dx = LogDensityProblems.logdensity_and_gradient(tpost, x)
        g .= -dx
        return nothing
    end

    before = objective(x0)
    initial_loglike = Comrade.loglikelihood(post, initial)

    # Optimize in the transformed coordinates; inspect convergence below.
    println("basis=sym6; setup_seconds=", time() - started)
    optstart = time()
    result = Optim.optimize(objective, gradient!, x0, Optim.LBFGS(), Optim.Options(; iterations = maxiters, show_trace = true))
    fitted = Comrade.transform(tpost, Optim.minimizer(result))
    after = objective(Optim.minimizer(result))
    fitted_model = skymodel(post, fitted)
    fitted_image = intensitymap(fitted_model, grid)

    # Save the full model and the compact raster on the inference grid.
    save_fits(joinpath(output, "wavelet_closure_map.fits"), fitted_image)
    fitted_fluctuation = apply(synthesis, reshape(fitted.sky.l, dims))
    fitted_raster = (1 - fitted.sky.fg) .* apply_fluctuations(CenteredLR(), mean_image, fitted.sky.σimg .* fitted_fluctuation)
    save_fits(joinpath(output, "wavelet_closure_raster.fits"), fitted_raster)

    # Match Comrade's tutorial: evaluate the continuous fitted model on a
    # separate 100×100 display grid without changing the inference grid.
    display_grid = imagepixels(μas2rad(150.0), μas2rad(150.0), 100, 100)
    display_image = intensitymap(fitted_model, display_grid)
    save_fits(joinpath(output, "wavelet_closure_map_display100.fits"), display_image)

    # Report optimizer status and fluxes alongside the saved file paths.
    println("UVFITS: ", data_path, "; grid=", dims, "; wavelet levels=", levels, "; basis=sym6; seed=", seed)
    println("optimization_and_output_seconds=", time() - optstart, "; total_seconds=", time() - started)
    println("objective initial=", before, "; final=", after, "; iterations=", Optim.iterations(result), "; converged=", Optim.converged(result))
    println("data log likelihood initial=", initial_loglike, "; final=", Comrade.loglikelihood(post, fitted))
    println("iteration limit reached=", Optim.iteration_limit_reached(result))
    println(
        "image min=", minimum(baseimage(fitted_image)), "; sampled image flux=", Comrade.flux(fitted_image),
        "; model total flux=", Comrade.flux(fitted_model)
    )
    println("σimg=", fitted.sky.σimg, "; raster total flux=", 1 - fitted.sky.fg, "; extended total flux=", fitted.sky.fg)
    println("image saved: ", joinpath(output, "wavelet_closure_map.fits"))
    println("inference-grid raster saved: ", joinpath(output, "wavelet_closure_raster.fits"))
    println("100×100 display-grid model saved: ", joinpath(output, "wavelet_closure_map_display100.fits"))
    return (; post, fitted, fitted_image, result, before, after)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
