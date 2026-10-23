@testset "Independent latent distributions" begin
    for T in (Float32, Float64), dims in ((7,), (3, 4), (2, 3, 4), (2, 2, 2, 2))
        for d in (Normal(T(0.2), T(1.3)), Laplace(T(-0.1), T(0.7)))
            f = UnivariateRandomField(dims, d)
            @test size(f) == dims
            @test length(f) == prod(dims)
            @test ndims(f) == length(dims)
            @test size(f, length(dims) + 1) == 1
            @test eltype(f) == T
            @test eltype(typeof(f)) == T
            z = rand(MersenneTwister(23), f)
            @test eltype(z) == T
            @test insupport(f, z)
            @test !insupport(f, fill(T(NaN), dims))
            @test logpdf(f, z) ≈ sum(logpdf.(d, z))
            @test @inferred(logpdf(f, z)) isa Real
            a, b = similar(z), similar(z)
            @test rand!(MersenneTwister(12), f, a) === a
            rand!(MersenneTwister(12), d, b)
            @test a == b
            batch = Array{T}(undef, dims..., 3)
            rand!(MersenneTwister(1), f, batch)
            @test length(logpdf(f, batch)) == 3
            samples = rand(MersenneTwister(1), f, 3)
            @test ndims(f) == 1 ? size(samples) == (length(f), 3) : length(samples) == 3
            @test_throws DimensionMismatch logpdf(f, zeros(T, ntuple(_ -> 1, length(dims))))
        end
    end
    f = UnivariateGaussianRandomField((3,); μ=0.4, σ=1.2)
    z = [0.2, -0.7, 1.3]
    @test logpdf(f, z) ≈ logpdf(MvNormal(fill(0.4, 3), 1.2^2 * I), z)
    @test UnivariateLaplaceRandomField(Float32, (2,)).dist == Laplace(0.0f0, 1.0f0)
    @test UnivariateGaussianRandomField((2,), Normal()) ==
        UnivariateRandomField((2,), Normal())
    @test_throws ArgumentError UnivariateRandomField((0, 2), Normal())
    @test_throws ArgumentError UnivariateRandomField((), Normal())
    @test_throws DomainError UnivariateLaplaceRandomField((2,); b=0)
    @test_throws ArgumentError size(f, 0)
end
