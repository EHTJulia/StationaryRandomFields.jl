using Reactant, Enzyme, Distributions, StationaryRandomFields, Statistics

function transform_loss(x, p)
    return sum(abs2, apply(p, x) .- 0.3)
end
function transform_gradient(x, p)
    dx = zero(x)
    Enzyme.autodiff(Reverse, transform_loss, Active, Duplicated(x, dx), Const(p))
    return dx
end
function measure(label, f, args...)
    Reactant.synchronize(f(args...))
    times = Float64[]
    bytes = Int[]
    for _ in 1:30
        trial = @timed Reactant.synchronize(f(args...))
        push!(times, trial.time)
        push!(bytes, trial.bytes)
    end
    println(
        label, " | median_us=", median(times) * 1e6, " | median_host_bytes=", median(bytes)
    )
    return nothing
end
function main()
    backend = get(ENV, "SRF_REACTANT_BACKEND", "cpu")
    Reactant.set_default_backend(backend)
    println(
        "Julia=",
        VERSION,
        " Reactant=",
        pkgversion(Reactant),
        " Enzyme=",
        pkgversion(Enzyme),
        " backend=",
        backend,
        " devices=",
        Reactant.devices(),
        " Julia threads=",
        Threads.nthreads(),
    )
    println(
        "samples=30, evals=1; resident inputs; synchronized; host bytes exclude device allocation",
    )
    for T in (Float32, Float64), dims in ((32, 32), (8, 8, 8))
        x = Reactant.to_rarray(reshape(T.(sin.(1:prod(dims))), dims))
        println("case=", T, " ", dims)
        for t in (
            HartleyTransform(dims),
            WaveletTransform(dims; levels=2),
            PeriodicCovariance(
                FieldGrid(dims), SinglePowerLaw{length(dims)}(T(-0.8); dc=T(0.3))
            ),
        )
            println("plan construction ", typeof(t))
            @time p = prepare(t, x)
            println("compile forward ", typeof(t))
            @time f = Reactant.@compile apply(p, x)
            println("compile gradient ", typeof(t))
            @time df = Reactant.@compile transform_gradient(x, p)
            measure("forward $(typeof(t))", f, p, x)
            measure("gradient $(typeof(t))", df, x, p)
        end
        latent = UnivariateGaussianRandomField(T, dims)
        @time density = Reactant.@compile logpdf(latent, x)
        measure("latent logpdf", density, latent, x)
        s = prepare(
            SpectralScaling(FieldGrid(dims), ConstantPowerSpectrum{length(dims)}(T(2))), x
        )
        @time scaling = Reactant.@compile apply(s, x)
        measure("spectrum apply", scaling, s, x)
    end
end
main()
