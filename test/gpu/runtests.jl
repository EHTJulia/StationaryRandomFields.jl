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
            f = UnivariateRandomField(dims, dist)
            @test logpdf(f, dx) ≈ logpdf(f, x)
            rand!(f, dx)
            @test all(isfinite, Array(dx))
        end
    end
    for dims in ((8,), (4, 8), (4, 4, 4), (4, 2, 2, 2))
        x = reshape(Float32.(sin.(1:prod(dims))), dims)
        for wt in (wavelet(WT.haar), wavelet(WT.db2))
            t = WaveletTransform(dims, wt; levels=1)
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
