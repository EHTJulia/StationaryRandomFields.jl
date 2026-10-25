@testset "Deferred transform preparation" begin
    for dims in ((16,), (16, 16)), T in (Float32, Float64)
        latent = UnivariateGaussianRandomField(T, dims)
        z = rand(MersenneTwister(14), latent)
        grid = FieldGrid(dims; spacing = ntuple(_ -> 0.2, length(dims)))
        model = MaternPS{length(dims)}(T(2), T(1.5))
        spec = TransformChain(SpectralScaling(model), HartleyTransform())
        p = prepare(spec, latent; grid)
        old = prepare(TransformChain(SpectralScaling(grid, model), HartleyTransform(dims)), z)
        @test apply(p, z) ≈ apply(old, z)
        @test apply(prepare(spec, z; grid), z) ≈ apply(p, z)
        @test apply(prepare(spec, latent; prototype = z, grid), z) ≈ apply(p, z)
        @test (@inferred apply(p, z)) ≈ apply(old, z)
        @test !applicable(logpdf, p, z)
        @test_throws ArgumentError prepare(spec, latent; grid, spacing = grid.spacing)
        @test_throws DimensionMismatch prepare(spec, latent; prototype = zeros(T, 3))
        wt = WaveletTransform(wavelet(WT.db4))
        @test apply(prepare(wt, latent), z) ≈ apply(prepare(WaveletTransform(dims, wavelet(WT.db4)), z), z)
        @test prepare(wt, latent).levels == 4
        normalized = prepare(SpectralScaling(model; normalization = :unit_variance), latent; grid)
        amplitudes = apply(normalized, ones(T, dims))
        @test sum(abs2, amplitudes) / length(amplitudes) ≈ one(T)
        changed = MaternPS{length(dims)}(T(3), T(2))
        @test sum(abs2, apply(normalized, ones(T, dims), changed)) / length(z) ≈ one(T)
        @test occursin("maximum", repr(wt))
        @test occursin("rho=", repr(spec))
    end
end

function normalized_spectrum_loss(rho, z, p)
    model = MaternPS{1}(rho, 1.5)
    return sum(abs2, apply(p, z, model) .- 0.3)
end
@testset "Normalized spectrum differentiation" begin
    latent = UnivariateGaussianRandomField((16,))
    z = rand(MersenneTwister(2), latent)
    p = prepare(SpectralScaling(MaternPS{1}(2.0, 1.5); normalization = :unit_variance), latent)
    rho = 2.0
    derivative = Enzyme.autodiff(Reverse, normalized_spectrum_loss, Active, Active(rho), Const(z), Const(p))[1][1]
    h = 1.0e-5
    reference = (normalized_spectrum_loss(rho + h, z, p) - normalized_spectrum_loss(rho - h, z, p)) / (2h)
    @test derivative ≈ reference rtol = 1.0e-5
end
