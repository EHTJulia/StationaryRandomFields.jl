using CUDA, Distributions, Random, Test, Wavelets, StationaryRandomFields

CUDA.functional(true) ||
    error("A working CUDA GPU is required; this suite must not silently skip.")
CUDA.allowscalar(false)
CUDA.versioninfo()

@testset "CUDA field operations" begin
    for T in (Float32, Float64), dims in ((7,), (3, 4), (2, 3, 4))
        x = reshape(T.(sin.(1:prod(dims))), dims)
        dx = CuArray(x)
        for t in (
            HartleyTransform(dims),
            SpectralScaling(FieldGrid(dims), ConstantPowerSpectrum{length(dims)}(T(2))),
        )
            p = prepare(t, dx)
            @test Array(apply(p, dx)) ≈ apply(prepare(t, x), x)
            y = similar(dx)
            @test apply!(y, p, dx) === y
            @test Array(y) ≈ apply(prepare(t, x), x)
        end
        for dist in (Normal(T(0.1), T(1.3)), Laplace(T(-0.2), T(0.7)))
            copyto!(dx, x)
            f = UnivariateRandomField(dims, dist)
            @test logpdf(f, dx) ≈ logpdf(f, x)
            rand!(f, dx)
            @test all(isfinite, Array(dx))
        end
    end
    for dims in ((8,), (4, 8), (4, 4, 4), (4, 2, 2, 2))
        x = reshape(Float32.(sin.(1:prod(dims))), dims)
        for wt in (wavelet(WT.haar), wavelet(WT.db2), wavelet(WT.db4))
            t = WaveletTransform(dims, wt)
            @test Array(apply(prepare(t, CuArray(x)), CuArray(x))) ≈ apply(prepare(t, x), x)
        end
    end
    f = UnivariateGaussianRandomField(Float32, (8, 8))
    field = TransformedRandomField(
        f,
        TransformChain(
            SpectralScaling(FieldGrid(size(f)), ConstantPowerSpectrum{2}(1.0f0)),
            HartleyTransform(size(f)),
        ),
    )
    output = CUDA.zeros(Float32, size(f))
    @test rand!(field, output) === output
    @test all(isfinite, Array(output))
end

@testset "CUDA covariance action" begin
    for T in (Float32, Float64), dims in ((7,), (3, 4), (2, 3, 4))
        x = reshape(T.(sin.(1:prod(dims))), dims)
        grid = FieldGrid(dims)
        model = SinglePowerLaw{length(dims)}(T(-0.8); dc=T(0.3))
        covariance = PeriodicCovariance(grid, model; variance=T(1.4))
        dx = CuArray(x)
        @test Array(apply(prepare(covariance, dx), dx)) ≈ apply(prepare(covariance, x), x)
    end
end

@testset "CUDA multilevel wavelet synthesis and transpose" begin
    for T in (Float32, Float64), dims in ((16,), (8, 16), (4, 8, 4), (4, 2, 2, 2), (24, 40), (7, 8), (1,))
        x = reshape(T.(sin.(1:prod(dims))), dims)
        v = reshape(T.(cos.(1:prod(dims))), dims)
        for wt in (wavelet(WT.haar), wavelet(WT.db4)), levels in unique((0, minimum(trailing_zeros, dims)))
            descriptor = WaveletTransform(dims, wt; levels)
            dx, dv = CuArray(x), CuArray(v)
            p, pc = prepare(descriptor, dx), prepare(descriptor, x)
            y = @inferred apply(p, dx)
            transpose_v = @inferred StationaryRandomFields.wavelet_analysis(p, dv)
            @test Array(y) ≈ apply(pc, x)
            @test Array(transpose_v) ≈ StationaryRandomFields.wavelet_analysis(pc, v)
            @test Array(StationaryRandomFields.wavelet_analysis(p, y)) ≈ x
            @test sum(Array(y) .* v) ≈ sum(x .* Array(transpose_v))
            @test Array(dx) == x
            @test Array(dv) == v
            output = similar(dx)
            @test apply!(output, p, dx) === output
            @test Array(output) ≈ Array(y)
            @test_throws DimensionMismatch apply(p, CUDA.zeros(T, prod(dims) + 1))
        end
    end
end

include("../backend_descriptors.jl")
check_backend_descriptors(CuArray, apply)
