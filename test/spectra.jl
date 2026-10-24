@testset "Power spectra" begin
    p = SaturatedPowerLaw{2}(-2.0f0, 0.1f0, 2.0f0)
    q = ConstantPowerSpectrum{2}(0.3f0)
    @test power_point(p, 0.0f0, 0.0f0) ≈ 100.0f0
    @test power_point(p, 3.0f0, 0.0f0) == 0.0f0
    @test power_point(p + q, 0.4f0, 0.2f0) ≈ power_point(p, 0.4f0, 0.2f0) + 0.3f0
    @test amplitude_point(4p, 0.4f0, 0.2f0) ≈ 2amplitude_point(p, 0.4f0, 0.2f0)
    @test power_point(p / 2, 0.4f0, 0.2f0) ≈ power_point(p, 0.4f0, 0.2f0) / 2
    @test power_point(stretched(p, 2, 1), 0.8, 0.2) ≈ power_point(p, 0.4, 0.2)
    θ = 0.37
    r = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    anisotropic = stretched(p, 3, 1)
    v = [0.3, 0.7]
    @test power_point(rotated(anisotropic, r), v...) ≈ power_point(anisotropic, (r * v)...)
    @test power_point(rotated(anisotropic, θ), v...) ≈
        power_point(rotated(anisotropic, r), v...)
    g = FieldGrid((5, 4); spacing = (0.2, 0.3))
    out = zeros(Float32, size(g))
    coords = SRF.frequency_coordinates(g, out)
    @test power_map!(out, p, coords) === out
    @test out ≈ power_map(p, g)
    @test amplitude_map!(out, p, coords) ≈ amplitude_map(p, g)
    @test @inferred(power_point(p, 0.3f0, 0.2f0)) isa Float32
    @test power_point(SinglePowerLaw{2}(-2), 0.0, 0.0) == 0
    @test_throws ArgumentError ConstantPowerSpectrum{2}(-1)
    @test_throws ArgumentError -2p
    @test_throws ArgumentError p / 0
    @test_throws ArgumentError SaturatedPowerLaw{2}(-2, 1, 0.1)
    @test_throws ArgumentError FieldGrid((2,); spacing = (-1.0,))
    @test_throws ArgumentError Rotate([1.0 1.0; 0.0 1.0])
    @test_throws DimensionMismatch rotated(p, Matrix{Float64}(I, 3, 3))
    @test_throws DimensionMismatch power_point(p, 1.0)
end

@testset "Matérn and Markov spectral shapes" begin
    for T in (Float32, Float64), N in (1, 2, 3)
        m = MaternPS{N}(T(2), T(1.5))
        p = MarkovPS{N}((T(2), T(0.5)))
        f = ntuple(_ -> T(0.1), N)
        k2 = (2pi)^2 * sum(abs2, f)
        @test power_point(m, f...) ≈ (1 + 4k2 / 12)^(-1.5 - N / 2)
        @test power_point(p, f...) ≈ inv(1 + 4k2 + (0.25k2)^2)
        for model in (m, p)
            @test power_point(model, ntuple(_ -> 0, N)...) == 1
            @test @inferred(power_point(model, f...)) isa T
            @test power_point(model, f...) == power_point(model, (-1 .* f)...)
            @test amplitude_point(model, f...)^2 ≈ power_point(model, f...)
            @test power_point(model + ConstantPowerSpectrum{N}(0.5), f...) ≈ power_point(model, f...) + 0.5 rtol = sqrt(eps(T))
            g = FieldGrid(ntuple(_ -> 4, N))
            x = ones(T, size(g))
            prepared = prepare(SpectralScaling(g, model), x)
            @test apply(prepared, x) ≈ amplitude_map(model, g)
            @test power_point(stretched(model, ntuple(_ -> 2, N)...), (2 .* f)...) ≈ power_point(model, f...)
        end
    end
    @test power_point(MaternPS{1}(2, 0.5), 0.2) ≈ power_point(MarkovPS{1}((1,)), 0.2)
    @test_throws ArgumentError MaternPS{2}(0, 1)
    @test_throws ArgumentError MaternPS{2}(1, 0)
    @test_throws ArgumentError MaternPS{0}(1, 1)
    @test_throws ArgumentError MarkovPS{2}(())
    @test_throws ArgumentError MarkovPS{2}((0, 0))
    @test_throws ArgumentError MarkovPS{2}((-1,))
    @test_throws ArgumentError MarkovPS{2}((Inf,))
    @test_throws DimensionMismatch power_point(MaternPS{2}(1, 1), 0.1)
end
