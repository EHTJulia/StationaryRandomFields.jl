# initialize some parameters
normal = Normal()
laplace = Laplace()

# check if initiralization works ok
@test UnivariateGaussianRandomField((10,), normal) == UnivariateRandomField((10,), normal)
@test UnivariateLaplaceRandomField((10,), laplace) == UnivariateRandomField((10,), laplace)

# check random numbers
grf = UnivariateGaussianRandomField((100, 100), normal)
lrf = UnivariateLaplaceRandomField((100, 100), laplace)

# check random number generations and logpdf
for rf in [grf, lrf]
    # check array generation
    @test zeros(rf) == zeros(eltype(rf), size(rf)...)
    @test ones(rf) == ones(eltype(rf), size(rf)...)
    @test fill(1, rf) == fill(eltype(rf)(1), size(rf)...)

    # check rand
    testarr1 = zeros(rf)
    testarr2 = zeros(rf)
    rand!(StableRNG(123), rf, testarr1)
    rand!(StableRNG(123), rf.dist, testarr2)
    @test testarr1 == testarr2

    # check logpdf
    @test logpdf(rf, testarr1) == sum(logpdf.(rf.dist, testarr2))
end