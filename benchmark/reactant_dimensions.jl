using Reactant, BenchmarkTools, StationaryRandomFields, Random, LinearAlgebra
using InteractiveUtils, Wavelets, Statistics

function measure(backend, operation, transform, cpu)
    x = Reactant.to_rarray(cpu)
    prepare_seconds = @elapsed p = prepare(transform, UnivariateGaussianRandomField(eltype(x), size(x)); prototype = x)
    compile_seconds = @elapsed compiled = Reactant.@compile apply(p, x)
    first_execution = @elapsed begin
        actual = compiled(p, x)
        Reactant.synchronize(actual)
    end
    expected = apply(prepare(transform, UnivariateGaussianRandomField(eltype(cpu), size(cpu)); prototype = cpu), cpu)
    tolerance = eltype(cpu) === Float32 ? 1.0e-4 : 1.0e-10
    isapprox(Array(actual), expected; rtol = tolerance) || error("CPU/Reactant mismatch")
    Array(x) == cpu || error("Input was modified")
    println("inferred_return,", backend, ',', operation, ',', eltype(cpu), ',', join(size(cpu), "x"), ',', Base.return_types(compiled, Tuple{typeof(p), typeof(x)}))
    for _ in 1:5
        Reactant.synchronize(compiled(p, x))
    end
    trial = @benchmark Reactant.synchronize($compiled($p, $x)) samples = 100 evals = 1 seconds = 2
    result = median(trial)
    Array(x) == cpu || error("Benchmark modified input")
    println(join((backend, operation, eltype(cpu), join(size(cpu), "x"), result.time / 1.0e6, result.memory, result.allocs, prepare_seconds, compile_seconds, first_execution, std(trial.times) / 1.0e6), ','))
    flush(stdout)
    return nothing
end

function main()
    backend = get(ENV, "SRF_REACTANT_BACKEND", "cpu")
    Reactant.set_default_backend(backend)
    BLAS.set_num_threads(1)
    fftw = Base.loaded_modules[Base.PkgId(Base.UUID("7a1cc6ca-52ef-59f5-83cd-3a7055c09341"), "FFTW")]
    fftw.set_num_threads(4)
    versioninfo()
    println("CPU=", Sys.cpu_info()[1].model)
    println("Reactant=", pkgversion(Reactant), " Wavelets=", pkgversion(Wavelets), " BenchmarkTools=", pkgversion(BenchmarkTools))
    println("SRF=", pkgversion(StationaryRandomFields), " FHT=", pkgversion(StationaryRandomFields.FastHartleyTransform))
    println("backend=", backend, " devices=", Reactant.devices(), " Julia_threads=", Threads.nthreads(), " FFTW_threads=", fftw.get_num_threads(), " BLAS_threads=", BLAS.get_num_threads())
    println("Reactant default compiler settings; resident arrays; unchanged input; synchronized; 5 warmups; 100 samples, evals=1, seconds=2; bytes are host allocations")
    println("backend,operation,type,shape,median_ms,host_bytes,host_allocs,prepare_s,compile_s,first_execution_s,std_ms")
    flush(stdout)
    for T in (Float32, Float64), dims in ((4096,), (256, 256), (64, 64, 64))
        cpu = randn(MersenneTwister(12), T, dims)
        spectrum = SinglePowerLaw{length(dims)}(T(-2); dc = T(0))
        transforms = (
            ("FHT", HartleyTransform()),
            ("Spectral synthesis", TransformChain(SpectralScaling(spectrum), HartleyTransform())),
            ("Db4 synthesis", WaveletTransform(wavelet(WT.db4))),
        )
        for (name, transform) in transforms
            println("starting,", backend, ',', name, ',', T, ',', dims)
            flush(stdout)
            measure(backend, name, transform, cpu)
        end
    end
    return nothing
end
main()
