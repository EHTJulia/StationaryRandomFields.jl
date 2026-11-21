@testset "Power spectra" begin
    for T in (Float32, Float64), n in (4, 5)
        axis = only(SRF.frequency_axes(FieldGrid((n,); spacing = (0.25,)), T))
        expected = n == 4 ? (0, 1, -2, -1) : (0, 1, 2, -2, -1)
        @test eltype(axis) == T
        @test axis ≈ T.(collect(expected)) .* (T(4) / n)
    end
    p = SaturatedPowerLaw{2}(-2.0f0, 0.1f0, 2.0f0)
    q = ConstantPowerSpectrum{2}(0.3f0)
    @test power_point(p, 0.0f0, 0.0f0) ≈ 100.0f0
    @test power_point(p, 3.0f0, 0.0f0) == 0.0f0
    @test power_point(p + q, 0.4f0, 0.2f0) ≈ power_point(p, 0.4f0, 0.2f0) + 0.3f0
    @test amplitude_point(4p, 0.4f0, 0.2f0) ≈ 2amplitude_point(p, 0.4f0, 0.2f0)
    @test power_point(p / 2, 0.4f0, 0.2f0) ≈ power_point(p, 0.4f0, 0.2f0) / 2
    @test power_point(stretched_ps(p, 2, 1), 0.8, 0.2) ≈ power_point(p, 0.4, 0.2)
    θ = 0.37
    r = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    anisotropic = stretched_ps(p, 3, 1)
    v = [0.3, 0.7]
    @test power_point(rotated_ps(anisotropic, r), v...) ≈ power_point(anisotropic, (r * v)...)
    @test power_point(rotated_ps(anisotropic, θ), v...) ≈
        power_point(rotated_ps(anisotropic, r), v...)
    g = FieldGrid((5, 4); spacing = (0.2, 0.3))
    out = zeros(Float32, size(g))
    coords = SRF.frequency_coordinates(g, out)
    @test power_map!(out, p, coords) === out
    @test out ≈ power_map(p, g)
    @test amplitude_map!(out, p, coords) ≈ amplitude_map(p, g)
    @test @inferred(power_point(p, 0.3f0, 0.2f0)) isa Float32
    @test power_point(SinglePowerLaw{2}(-2), 0.0, 0.0) == 0

    @testset "SinglePowerLaw frequency precision" begin
        for T in (Float32, Float64)
            model = SinglePowerLaw{2}(T(-2))
            @test power_point(model, 1, 2) isa T
            @test power_point(model, 0, 0) isa T
            @test amplitude_point(model, 1, 2) isa T
            @test power_point(model, 1.0f0, 2.0f0) isa promote_type(T, Float32)
            @test power_point(model, 1.0, 2.0) isa promote_type(T, Float64)
        end
    end
    @test_throws ArgumentError ConstantPowerSpectrum{2}(-1)
    @test_throws ArgumentError -2p
    @test_throws ArgumentError p / 0
    @test_throws ArgumentError SaturatedPowerLaw{2}(-2, 1, 0.1)
    @test_throws ArgumentError FieldGrid((2,); spacing = (-1.0,))
    @test PowerSpectrumRotate([1.0 1.0; 0.0 1.0]).matrix == [1.0 1.0; 0.0 1.0]
    p3 = stretched_ps(SinglePowerLaw{3}(-2.0), 1.0, 2.0, 3.0)
    r3 = [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0]
    v3 = [0.2, 0.3, 0.4]
    @test @inferred(power_point(rotated_ps(p3, r3), v3...)) ≈ power_point(p3, (r3 * v3)...)
    @test_throws DimensionMismatch PowerSpectrumRotate(zeros(2, 3))
    @test_throws DimensionMismatch PowerSpectrumRotate(zeros(0, 0))
    @test_throws DimensionMismatch rotated_ps(p, Matrix{Float64}(I, 3, 3))
    modifiers = (PowerSpectrumRenormalize(2.0), PowerSpectrumStretch(3.0, 1.0), PowerSpectrumRotate(θ))
    mixed = @inferred modify_ps(p, modifiers...)
    @test power_point(mixed, v...) ≈ power_point(rotated_ps(stretched_ps(renormed_ps(p, 2.0), 3.0, 1.0), θ), v...)
    @test @inferred(modify_ps(p)) isa SRF.ModifiedPowerSpectrumModel
    @test_throws DimensionMismatch modify_ps(p, modifiers[1], PowerSpectrumStretch(1.0, 1.0, 1.0))
    @test_throws DimensionMismatch modify_ps(p, modifiers[1], PowerSpectrumRotate(Matrix{Float64}(I, 3, 3)))
    @test_throws MethodError power_point(p, 1.0)
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
            @test power_point(stretched_ps(model, ntuple(_ -> 2, N)...), (2 .* f)...) ≈ power_point(model, f...)
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
    @test_throws MethodError power_point(MaternPS{2}(1, 1), 0.1)
end

@testset "Spectrum frequency dimensions" begin
    models = (
        ConstantPowerSpectrum{2}(1.0),
        SinglePowerLaw{2}(-2.0),
        SaturatedPowerLaw{2}(-2.0, 0.1, 2.0),
        MaternPS{2}(1.0, 1.0),
        MarkovPS{2}((1.0,)),
    )
    for model in (models..., models[1] + models[2], stretched_ps(models[1], 1.0, 2.0))
        @test power_point(model, 0.2, 0.3) isa Real
        @test amplitude_point(model, 0.2, 0.3) isa Real
        for coordinates in ((), (0.2,), (0.2, 0.3, 0.4))
            @test_throws MethodError power_point(model, coordinates...)
            @test_throws MethodError amplitude_point(model, coordinates...)
        end
        @test_throws DimensionMismatch SRF.paired_power(model, 0.2, 0.3)
        @test_throws DimensionMismatch SRF.paired_power(model, 0.2, 0.3, -0.2, -0.3, 0.4)
    end
end
