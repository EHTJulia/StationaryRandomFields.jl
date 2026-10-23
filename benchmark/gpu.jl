using CUDA, BenchmarkTools, StationaryRandomFields

function main()
    CUDA.functional(true) || error("A working CUDA device is required.")
    CUDA.allowscalar(false)
    CUDA.versioninfo()
    println(
        "SRF=",
        pkgversion(StationaryRandomFields),
        " FHT=",
        pkgversion(StationaryRandomFields.FastHartleyTransform),
        " Wavelets=",
        pkgversion(StationaryRandomFields.Wavelets),
        " Julia threads=",
        Threads.nthreads(),
    )
    println(
        "Samples=30, evals=1; each execution synchronizes; execution excludes input/output transfers; preparation includes metadata copies; reported bytes are host allocations.",
    )
    for T in (Float32, Float64), dims in ((256, 256), (64, 64, 64))
        x = CUDA.randn(T, dims)
        println("case=", T, " ", dims)
        for t in (
            HartleyTransform(dims),
            WaveletTransform(dims; levels=2),
            PeriodicCovariance(
                FieldGrid(dims), SinglePowerLaw{length(dims)}(T(-0.8); dc=T(0.3))
            ),
        )
            println("first preparation (includes compilation) ", typeof(t))
            @time CUDA.@sync prepare(t, x)
            println("plan ", typeof(t), " ", dims)
            display(@benchmark CUDA.@sync(prepare($t, $x)) samples=30 evals=1 seconds=1)
            p = prepare(t, x)
            println("first execution (includes compilation) ", typeof(t))
            @time CUDA.@sync apply(p, x)
            println("execution ", typeof(t), " ", dims)
            display(@benchmark CUDA.@sync(apply($p, $x)) samples=30 evals=1 seconds=1)
        end
    end
end
main()
