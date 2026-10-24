using CUDA, Wavelets, StationaryRandomFields, BenchmarkTools, Random
function main()
    CUDA.allowscalar(false)
    CUDA.versioninfo()
    println("Wavelets GPU extension: ", Base.get_extension(Wavelets, :WaveletsGPUExt))
    println("type,shape,levels,path,ms,host_bytes,host_allocs")
    for T in (Float32, Float64), dims in ((4096,), (256, 256), (64, 64, 64))
        x = randn(MersenneTwister(12), T, dims)
        dx = CuArray(x)
        p = prepare(WaveletTransform(dims, wavelet(WT.db4)), dx)
        expected = Wavelets.idwt(x, p.wavelet, p.levels)
        for (name, f) in (("generic", () -> StationaryRandomFields.wavelet_synthesis(dx, p.filters, p.levels)), ("wavelets", () -> Wavelets.idwt(dx, p.wavelet, p.levels)))
            y = f()
            @assert isapprox(Array(y), expected; rtol = T === Float32 ? 1.0e-4 : 1.0e-10)
            @assert Array(dx) == x
            for _ in 1:5
                f(); CUDA.synchronize()
            end
            trial = @benchmark ($f(); CUDA.synchronize()) samples = 100 evals = 1 seconds = 2
            m = median(trial)
            println(join((T, join(dims, "x"), p.levels, name, m.time / 1.0e6, m.memory, m.allocs), ','))
        end
    end
    return
end
main()
