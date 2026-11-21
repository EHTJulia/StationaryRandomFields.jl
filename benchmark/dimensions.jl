using CUDA, BenchmarkTools, StationaryRandomFields, Random, LinearAlgebra
using InteractiveUtils, Wavelets, Statistics

function measure(backend, operation, transform, x)
    latent = UnivariateGaussianRandomField(eltype(x), size(x))
    sync() = backend == "CUDA" ? CUDA.synchronize() : nothing
    sync()
    first_prepare = @elapsed p = prepare(transform, UnivariateGaussianRandomField(eltype(x), size(x)); prototype = x)
    sync()
    p = prepare(transform, UnivariateGaussianRandomField(eltype(x), size(x)); prototype = x)
    sync()
    planning = @benchmark (prepare($transform, $latent; prototype = $x); $sync()) samples = 30 evals = 1 seconds = 1
    first_execution = @elapsed begin
        apply(p, x)
        sync()
    end
    for _ in 1:5
        apply(p, x)
        sync()
    end
    trial = @benchmark (apply($p, $x); $sync()) samples = 100 evals = 1 seconds = 2
    result = median(trial)
    println(
        join(
            (
                backend, operation, eltype(x), join(size(x), "x"),
                result.time / 1.0e6, result.memory, result.allocs,
                median(planning).time / 1.0e6, first_prepare, first_execution, std(trial.times) / 1.0e6,
            ), ','
        )
    )
    return nothing
end

function main()
    CUDA.functional(true) || error("CUDA GPU required")
    CUDA.allowscalar(false)
    versioninfo()
    CUDA.versioninfo()
    BLAS.set_num_threads(1)
    fftw = Base.loaded_modules[Base.PkgId(Base.UUID("7a1cc6ca-52ef-59f5-83cd-3a7055c09341"), "FFTW")]
    fftw.set_num_threads(4)
    println("CPU=", Sys.cpu_info()[1].model)
    println("SRF=", pkgversion(StationaryRandomFields), " FHT=", pkgversion(StationaryRandomFields.FastHartleyTransform), " FFTW=", pkgversion(fftw), " CUDA=", pkgversion(CUDA))
    println("Julia_threads=", Threads.nthreads(), " FFTW_threads=", fftw.get_num_threads(), " BLAS_threads=", BLAS.get_num_threads())
    println("ESTIMATE; resident arrays; unchanged input; synchronized CUDA; 5 warmups; 100 samples, evals=1, seconds=2; bytes are host allocations")
    println("backend,operation,type,shape,median_ms,host_bytes,host_allocs,prepare_ms,first_prepare_s,first_execution_s,std_ms")
    for T in (Float32, Float64), dims in ((4096,), (256, 256), (64, 64, 64))
        cpu = randn(MersenneTwister(12), T, dims)
        gpu = CuArray(cpu)
        spectrum = SinglePowerLaw{length(dims)}(T(-2); dc = T(0))
        transforms = (
            ("FHT", HartleyTransform()),
            ("Spectral synthesis", TransformChain(SpectralScaling(spectrum), HartleyTransform())),
            ("Db4 synthesis", WaveletTransform(wavelet(WT.db4))),
        )
        for (name, transform) in transforms
            pc = prepare(transform, UnivariateGaussianRandomField(eltype(cpu), size(cpu)); prototype = cpu)
            pg = prepare(transform, UnivariateGaussianRandomField(eltype(gpu), size(gpu)); prototype = gpu)
            expected = apply(pc, cpu)
            actual = Array(apply(pg, gpu))
            isapprox(actual, expected; rtol = T == Float32 ? 1.0e-4 : 1.0e-10) || error("CPU/GPU mismatch")
            measure("CPU", name, transform, cpu)
            measure("CUDA", name, transform, gpu)
        end
    end
    return
end
main()
