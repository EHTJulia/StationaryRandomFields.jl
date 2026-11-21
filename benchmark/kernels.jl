using BenchmarkTools,
    Distributions,
    Enzyme,
    FFTW,
    InteractiveUtils,
    LinearAlgebra,
    Random,
    StationaryRandomFields
using KernelFunctions: kernelmatrix

function report(label, f, args...)
    f(args...)
    trial = @benchmark $f($args...) samples = 30 evals = 1 seconds = 0.5
    estimate = median(trial)
    return println(
        label,
        " | median_us=",
        estimate.time / 1000,
        " | bytes=",
        estimate.memory,
        " | allocs=",
        estimate.allocs,
    )
end
function covariance_gradient(a, p, x)
    model = SinglePowerLaw{ndims(x)}(a; dc = 0.3)
    return sum(abs2, apply(p, x, model))
end
function index_gradient(a, p, x)
    return Enzyme.autodiff(
        Reverse, covariance_gradient, Active, Active(a), Const(p), Const(x)
    )[1][1]
end
function main()
    versioninfo()
    println(
        "Julia=",
        VERSION,
        " SRF=",
        pkgversion(StationaryRandomFields),
        " KernelFunctions=",
        pkgversion(parentmodule(kernelmatrix)),
        " FFTW=",
        pkgversion(FFTW),
        " Enzyme=",
        pkgversion(Enzyme),
    )
    println(
        "Julia threads=",
        Threads.nthreads(),
        " FFTW threads=",
        FFTW.get_num_threads(),
        " BLAS threads=",
        BLAS.get_num_threads(),
    )
    println(
        "FFTW default ESTIMATE; samples=30 evals=1 seconds=0.5; CPU synchronized by return; no input mutation",
    )
    for T in (Float32, Float64), dims in ((32,), (8, 8), (4, 4, 4))
        grid = FieldGrid(dims; spacing = ntuple(_ -> one(T), length(dims)))
        model = SinglePowerLaw{length(dims)}(T(-0.8); dc = T(0.3))
        c = PeriodicCovariance(grid, model)
        x = reshape(T.(sin.(1:prod(dims))), dims)
        points = [
            ntuple(d -> T(i[d] - 1), length(dims)) for
                i in vec(collect(CartesianIndices(dims)))
        ]
        println("case=", T, " ", dims)
        println("first preparation (includes compilation)")
        @time p = prepare(c, x)
        println("first covariance action (includes compilation)")
        @time apply(p, x)
        println("first spectrum derivative (includes compilation)")
        @time index_gradient(T(-0.8), p, x)
        println(
            "inferred point=",
            Core.Compiler.return_type(
                covariance_value,
                Tuple{typeof(c), typeof(first(points)), typeof(last(points))},
            ),
        )
        println(
            "inferred action=", Core.Compiler.return_type(apply, Tuple{typeof(p), typeof(x)})
        )
        report("prepare", prepare, c, x)
        report("point", covariance_value, c, first(points), last(points))
        report("dense kernelmatrix", kernelmatrix, spectral_kernel(c), points)
        report("FHT action", apply, p, x)
        report("index derivative", index_gradient, T(-0.8), p, x)
    end
    return
end
main()
