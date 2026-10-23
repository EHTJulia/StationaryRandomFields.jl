using CUDA, BenchmarkTools, StationaryRandomFields

function main()
    CUDA.functional(true) || error("A working CUDA device is required.")
    CUDA.allowscalar(false)
    CUDA.versioninfo()
    println("Samples=30, evals=1; each execution synchronizes; transfers are excluded.")
    for dims in ((256, 256), (64, 64, 64))
        x = CUDA.randn(Float32, dims)
        for t in (HartleyTransform(dims), WaveletTransform(dims; levels=2))
            CUDA.@sync prepare(t, x)
            println("plan ", typeof(t), " ", dims)
            display(@benchmark CUDA.@sync(prepare($t, $x)) samples=30 evals=1 seconds=1)
            p = prepare(t, x)
            CUDA.@sync apply(p, x)
            println("execution ", typeof(t), " ", dims)
            display(@benchmark CUDA.@sync(apply($p, $x)) samples=30 evals=1 seconds=1)
        end
    end
end
main()
