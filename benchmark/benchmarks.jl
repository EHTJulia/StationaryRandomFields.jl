using BenchmarkTools, StationaryRandomFields, Distributions, Random, Enzyme, FFTW
using InteractiveUtils

function gradient!(dx, x, f)
    fill!(dx, zero(eltype(dx)))
    Enzyme.autodiff(Reverse, logpdf, Active, Const(f), Duplicated(x, dx))
    return dx
end
function report(name, b)
    trial = run(b; samples=30, evals=1, seconds=0.5)
    m = median(trial)
    println(name, " | ", round(m.time / 1000; digits=3), " | ", m.memory, " | ", m.allocs)
end
function main()
    versioninfo()
    println(
        "SRF=",
        pkgversion(StationaryRandomFields),
        " Distributions=",
        pkgversion(Distributions),
        " Enzyme=",
        pkgversion(Enzyme),
        " FFTW=",
        pkgversion(FFTW),
    )
    println(
        "FHT=",
        pkgversion(StationaryRandomFields.FastHartleyTransform),
        " Wavelets=",
        pkgversion(StationaryRandomFields.Wavelets),
        " BenchmarkTools=",
        pkgversion(BenchmarkTools),
    )
    println("FFTW threads=", FFTW.get_num_threads(), " Julia threads=", Threads.nthreads())
    println("Planner=ESTIMATE; samples=30; evals=1; seconds=0.5; median; CPU synchronous")
    println("operation | median μs | bytes | allocations")
    for T in (Float32, Float64), dims in ((64, 64), (16, 16, 16))
        println("\n", T, " ", dims)
        rng = MersenneTwister(12)
        f = UnivariateGaussianRandomField(T, dims)
        l = UnivariateLaplaceRandomField(T, dims)
        x = rand(rng, f)
        println(
            "First FHT call (planning + compilation + execution; @time reports compilation):",
        )
        @time apply(HartleyTransform(dims), x)
        y, dx = similar(x), zero(x)
        s = prepare(
            SpectralScaling(FieldGrid(dims), ConstantPowerSpectrum{length(dims)}(T(2))), x
        )
        h = prepare(HartleyTransform(dims), x)
        w = prepare(WaveletTransform(dims; levels=2), x)
        for (name, b) in (
            ("Gaussian sampling", @benchmarkable rand!($rng, $f, $x)),
            ("Laplace sampling", @benchmarkable rand!($rng, $l, $x)),
            ("Gaussian logpdf", @benchmarkable logpdf($f, $x)),
            ("Laplace logpdf", @benchmarkable logpdf($l, $x)),
            ("Spectrum apply", @benchmarkable apply($s, $x)),
            (
                "FHT plan construction",
                @benchmarkable prepare($(HartleyTransform(dims)), $x)
            ),
            ("FHT apply", @benchmarkable apply($h, $x)),
            ("FHT apply!", @benchmarkable apply!($y, $h, $x)),
            (
                "Wavelet plan construction",
                @benchmarkable prepare($(WaveletTransform(dims; levels=2)), $x)
            ),
            ("Wavelet apply", @benchmarkable apply($w, $x)),
            ("Latent gradient", @benchmarkable gradient!($dx, $x, $f)),
        )
            report(name, b)
        end
    end
end

main()
