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
    g = FieldGrid((5, 4); spacing=(0.2, 0.3))
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
    @test_throws ArgumentError FieldGrid((2,); spacing=(-1.0,))
    @test_throws ArgumentError Rotate([1.0 1.0; 0.0 1.0])
    @test_throws DimensionMismatch rotated(p, Matrix{Float64}(I, 3, 3))
    @test_throws DimensionMismatch power_point(p, 1.0)
end
