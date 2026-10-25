function direct_hartley(dims, ::Type{T}) where {T}
    inds = collect(CartesianIndices(dims))
    return (a -> reshape(a, prod(dims), prod(dims)))(
        [
            T(
                cospi(2sum((i[d] - 1) * (j[d] - 1) / dims[d] for d in eachindex(dims))) +
                    sinpi(2sum((i[d] - 1) * (j[d] - 1) / dims[d] for d in eachindex(dims))),
            ) for i in inds, j in inds
        ]
    )
end

@testset "Hartley synthesis and signal sampling" begin
    for T in (Float32, Float64), dims in ((7,), (3, 4), (2, 3, 2), (2, 2, 2, 2))
        z = reshape(T.(sin.(1:prod(dims))), dims)
        h = prepare(HartleyTransform(dims), z)
        q = direct_hartley(dims, T) / sqrt(T(length(z)))
        @test vec(apply(h, z)) ≈ q * vec(z)
        @test apply(h, apply(h, z)) ≈ z
        y = similar(z)
        @test apply!(y, h, z) === y
        @test y ≈ apply(h, z)
        @test @inferred(apply(h, z)) isa Array{T}
        @test_throws ArgumentError apply!(z, h, z)
        grid = FieldGrid(dims)
        s = prepare(SpectralScaling(grid, ConstantPowerSpectrum{length(dims)}(T(4))), z)
        @test apply(s, z) ≈ 2z
        @test apply(s, z, ConstantPowerSpectrum{length(dims)}(T(9))) ≈ 3z
        chain = TransformChain(s, h, ScaleTransform(dims, T(0.5)))
        @test apply(chain, z) ≈ apply(h, z)
        f = TransformedRandomField(UnivariateGaussianRandomField(T, dims), chain)
        @test size(rand(MersenneTwister(1), f)) == dims
        @test_throws ArgumentError logpdf(f, z)
        @test_throws ArgumentError insupport(f, z)
    end
    dims = (4, 5)
    z = zeros(dims)
    grid = FieldGrid(dims)
    p = rotated(stretched(SaturatedPowerLaw{2}(-1.3, 0.1, 1.0), 2.0, 1.0), 0.4)
    s = prepare(SpectralScaling(grid, p), z)
    h = prepare(HartleyTransform(dims), z)
    a = SRF.spectral_amplitude(s, p)
    q = direct_hartley(dims, Float64) / sqrt(length(z))
    expected = q * Diagonal(vec(a) .^ 2) * q'
    basis = Matrix{Float64}(I, length(z), length(z))
    generated = hcat(
        [vec(apply(h, apply(s, reshape(basis[:, i], dims)))) for i in 1:length(z)]...
    )
    @test generated * generated' ≈ expected
    indices = vec(collect(CartesianIndices(dims)))
    fourier = [
        cispi(-2sum((k[d] - 1) * (n[d] - 1) / dims[d] for d in eachindex(dims))) for
            k in indices, n in indices
    ]
    fourier_covariance = real.(fourier' * Diagonal(vec(a) .^ 2) * fourier) / length(z)
    @test expected ≈ fourier_covariance
    expected_periodogram = vec(sum(abs2, fourier * generated; dims = 2)) / length(z)
    @test expected_periodogram ≈ vec(a) .^ 2

    @test_throws DimensionMismatch TransformChain(
        HartleyTransform((3,)), HartleyTransform((4,))
    )
    @test_throws DimensionMismatch TransformedRandomField(
        UnivariateGaussianRandomField((3,)), HartleyTransform((4,))
    )
    @test_throws ArgumentError TransformChain()
end

@testset "Wavelet synthesis" begin
    for T in (Float32, Float64), dims in ((8,), (4, 8), (4, 4, 4)), levels in (0, 1, 2)
        x = reshape(T.(sin.(1:prod(dims))), dims)
        for wt in (wavelet(WT.haar), wavelet(WT.db2), wavelet(WT.coif2))
            t = prepare(WaveletTransform(dims, wt; levels), x)

            expected = idwt(x, wt, levels)
            actual = apply(t, x)

            @test actual ≈ expected
            @test SRF.wavelet_analysis(t, actual) ≈ x
            @test sum(abs2, actual) ≈ sum(abs2, x)
            @test @inferred(apply(t, x)) isa Array{T}
        end
    end

    x = reshape(sin.(1:256), 4, 4, 4, 4)
    t = prepare(WaveletTransform(size(x); levels = 2), x)
    @test SRF.wavelet_analysis(t, apply(t, x)) ≈ x
    @test_throws ArgumentError WaveletTransform((7, 8); levels = 1)
    @test_throws ArgumentError WaveletTransform((8,); levels = -1)
end

@testset "Maximum-depth periodic DWT" begin
    for T in (Float32, Float64), (dims, depth) in (
                ((16,), 4), ((8, 16), 3), ((8, 8, 8), 3),
                ((24, 40), 3), ((7, 8), 0), ((1,), 0),
            )

        x = reshape(T.(sin.(1:prod(dims))), dims)
        wt = wavelet(WT.db4)
        descriptor = @inferred WaveletTransform(dims, wt)
        @test descriptor.levels == depth
        p = prepare(descriptor, x)
        y = @inferred apply(p, x)
        @test y ≈ idwt(x, wt, depth)
        @test SRF.wavelet_synthesis(x, p.filters, depth) ≈ y
        @test SRF.wavelet_analysis(p, y) ≈ x
        @test sum(abs2, y) ≈ sum(abs2, x)
        @test WaveletTransform(dims, wt; levels = 0).levels == 0
    end
    @test WaveletTransform((16,)).levels == 4
    @test WaveletTransform((16,); levels = 1).levels == 1
    @test_throws ArgumentError WaveletTransform((8,); levels = 4)
end

struct ReshapeTransform{I, O} <: AbstractFieldTransform
    input::I
    output::O
end
StationaryRandomFields.input_size(t::ReshapeTransform) = t.input
StationaryRandomFields.output_size(t::ReshapeTransform) = t.output
StationaryRandomFields.prepare(t::ReshapeTransform, x) = t
StationaryRandomFields.apply(t::ReshapeTransform, x) = reshape(copy(x), t.output)

@testset "Shape-changing custom synthesis and composed sampling" begin
    latent = UnivariateLaplaceRandomField((12,))
    t = ReshapeTransform((12,), (3, 4))
    field = TransformedRandomField(latent, t)
    @test size(field) == (3, 4)
    @test eltype(typeof(field)) == Float64
    z = rand(MersenneTwister(5), latent)
    @test size(apply(t, z)) == (3, 4)
    @test rand(MersenneTwister(5), field) == reshape(z, 3, 4)
    @test length(rand(MersenneTwister(5), field, 3)) == 3
    batch = zeros(3, 4, 3)
    @test rand!(MersenneTwister(5), field, batch) === batch
    @test_throws ArgumentError logpdf(field, batch)
    grid = FieldGrid((4, 4))
    s = prepare(SpectralScaling(grid, SinglePowerLaw{2}(-2.0)), zeros(4, 4))
    @test @inferred(apply(s, ones(4, 4))) isa Matrix{Float64}
    @test apply(s, ones(4, 4))[1] == 0
end

@testset "Transform construction and element types" begin
    @test_throws ArgumentError IdentityTransform((0,))
    @test_throws ArgumentError HartleyTransform((0,))
    h = prepare(HartleyTransform((4,)), ones(Float32, 4))
    @test_throws ArgumentError apply(h, ones(Float64, 4))
    @test_throws DimensionMismatch ScaleTransform((3, 4), ones(4, 3))
    x = ones(Float32, 3, 4)
    @test eltype(apply(ScaleTransform(size(x), 1.5), x)) == Float32
    t = prepare(SpectralScaling(FieldGrid(size(x)), ConstantPowerSpectrum{2}(2.0)), x)
    @test eltype(apply(t, x)) == Float32
end
