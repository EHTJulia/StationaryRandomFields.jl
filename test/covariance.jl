using KernelFunctions: KernelFunctions
using KernelFunctions: kernelmatrix, kernelmatrix!, kernelmatrix_diag

function grid_points(grid)
    return [
        ntuple(d -> (i[d] - 1) * grid.spacing[d], ndims(grid)) for
        i in vec(collect(CartesianIndices(size(grid))))
    ]
end
function synthesis_matrix(grid, model, T)
    x = zeros(T, size(grid))
    transform = prepare(
        TransformChain(SpectralScaling(grid, model), HartleyTransform(size(grid))), x
    )
    return reduce(
        hcat,
        map(1:length(x)) do j
            basis = zero(x)
            basis[j] = 1
            vec(apply(transform, basis))
        end,
    )
end

@testset "Finite periodic covariance and KernelFunctions" begin
    for T in (Float32, Float64), dims in ((5,), (4,), (3, 4), (2, 3, 2), (2, 2, 2, 2))
        grid = FieldGrid(dims; spacing=ntuple(d -> T(0.7 + 0.1d), length(dims)))
        power =
            SinglePowerLaw{length(dims)}(T(-0.8); dc=T(0.3)) +
            ConstantPowerSpectrum{length(dims)}(T(0.2))
        latent = UnivariateGaussianRandomField(T, dims; σ=T(1.3))
        covariance = PeriodicCovariance(grid, power, latent)
        k = spectral_kernel(covariance)
        @test k isa KernelFunctions.Kernel
        points = grid_points(grid)
        matrix = kernelmatrix(k, points)
        synthesis = synthesis_matrix(grid, power, T)
        @test matrix ≈ var(latent.dist) .* (synthesis * synthesis') rtol=100eps(T)
        @test matrix ≈ matrix'
        @test eigmin(Symmetric(matrix)) >= -100eps(T) * opnorm(matrix)
        @test kernelmatrix_diag(k, points) ≈ diag(matrix)
        buffer = similar(matrix)
        @test kernelmatrix!(buffer, k, points) === buffer
        @test buffer ≈ matrix
        coordinates = reduce(hcat, collect.(points))
        @test kernelmatrix(k, coordinates; obsdim=2) ≈ matrix
        @test kernelmatrix(k, coordinates'; obsdim=1) ≈ matrix
        x = reshape(T.(sin.(1:prod(dims))), dims)
        p = prepare(covariance, x)
        @test vec(@inferred apply(p, x)) ≈ matrix * vec(x) rtol=100eps(T)
        y = similar(x)
        @test apply!(y, p, x) === y
        @test y ≈ apply(p, x)
        @test @inferred(covariance(first(points), last(points))) ≈ matrix[1, end]
        periods = size(grid) .* grid.spacing
        @test k(first(points) .+ periods, last(points)) ≈ matrix[1, end] atol=100eps(T)
        if length(dims) == 1
            @test k(first(points)[1], last(points)[1]) ≈ matrix[1, end]
        end
        @test_throws DimensionMismatch k(zeros(length(dims) + 1), zeros(length(dims)))
    end
    grid = FieldGrid((4, 5); spacing=(0.8, 1.1))
    model = rotated(stretched(SinglePowerLaw{2}(-1.2; dc=0.4), 1.8, 0.7), 0.4)
    c = PeriodicCovariance(grid, model)
    points = [(0.1sin(i), 0.7cos(2i)) for i in 1:9]
    matrix = kernelmatrix(spectral_kernel(c), points)
    @test matrix ≈ matrix'
    @test eigmin(Symmetric(matrix)) >= -1e-12 * opnorm(matrix)
    @test matrix ≈ [c(x, y) for x in points, y in points]
    ongrid = grid_points(grid)
    synthesis = synthesis_matrix(grid, model, Float64)
    @test kernelmatrix(spectral_kernel(c), ongrid) ≈ synthesis * synthesis'
    @test_throws DimensionMismatch PeriodicCovariance(grid, ConstantPowerSpectrum{1}(1))
    @test_throws ArgumentError PeriodicCovariance(grid, model; variance=-1)
    @test_throws DimensionMismatch PeriodicCovariance(
        grid, model, UnivariateGaussianRandomField((4, 4))
    )
    @test_throws DimensionMismatch c(0.1, 0.2)
    zero_covariance = PeriodicCovariance(grid, model; variance=0)
    @test iszero(zero_covariance((0.1, 0.2), (0.3, 0.4)))
end

@testset "Sampler empirical covariance" begin
    rng = MersenneTwister(283)
    grid = FieldGrid((4,))
    model = SinglePowerLaw{1}(-0.7; dc=0.5)
    latent = UnivariateGaussianRandomField((4,); σ=1.2)
    covariance = PeriodicCovariance(grid, model, latent)
    transform = TransformChain(SpectralScaling(grid, model), HartleyTransform(size(grid)))
    field = prepare(TransformedRandomField(latent, transform), zeros(4))
    draws = zeros(4, 20000)
    for j in axes(draws, 2)
        rand!(rng, field, view(draws, :, j))
    end
    empirical = draws * draws' / size(draws, 2)
    expected = kernelmatrix(spectral_kernel(covariance), grid_points(grid))
    @test norm(empirical - expected) / norm(expected) < 0.04
end

function covariance_model(a)
    return modify(
        SinglePowerLaw{2}(sum(a[2:2]); dc=0.3, check_args=false),
        Renormalize(exp(sum(a[1:1])); check_args=false),
        Stretch(exp(sum(a[3:3])), 0.8; check_args=false),
        Rotate(sum(a[4:4])),
    )
end
function kernel_parameter_loss(a, grid)
    c = PeriodicCovariance(grid, covariance_model(a))
    return spectral_kernel(c)((0.13, 0.21), (-0.27, 0.04))
end
function covariance_parameter_loss(a, p, x)
    return sum(abs2, apply(p, x, covariance_model(a)) .- 0.3)
end
@testset "Covariance spectrum parameter gradients" begin
    for T in (Float32, Float64)
        grid = FieldGrid((4, 5))
        a = T[0.1, -1.2, 0.2, 0.3]
        x = reshape(T.(sin.(1:20)), 4, 5)
        p = prepare(PeriodicCovariance(grid, covariance_model(a)), x)
        for (loss, args) in
            ((kernel_parameter_loss, (grid,)), (covariance_parameter_loss, (p, x)))
            da = zero(a)
            Enzyme.autodiff(Reverse, loss, Active, Duplicated(a, da), Const.(args)...)
            expected = finite_gradient(b -> loss(b, args...), a)
            @test all(isfinite, da)
            @test da ≈ expected rtol=(T == Float32 ? 0.003 : 1e-7)
            tangent = Enzyme.autodiff(
                Forward, loss, Duplicated, Duplicated(a, one.(a)), Const.(args)...
            )[1]
            @test tangent ≈ sum(expected) rtol=(T == Float32 ? 0.003 : 1e-7)
        end
    end
end

function mixed_index_loss(a, p, x)
    model = SinglePowerLaw{ndims(x)}(a; dc=0.3)
    return sum(abs2, apply(p, x, model))
end
@testset "Mixed precision covariance gradients" begin
    for T in (Float32, Float64), dims in ((8,), (3, 4))
        x = reshape(T.(sin.(1:prod(dims))), dims)
        grid = FieldGrid(dims)
        model = SinglePowerLaw{length(dims)}(T(-0.8); dc=T(0.3))
        p = prepare(PeriodicCovariance(grid, model), x)
        index = T(-0.8)
        actual = Enzyme.autodiff(
            Reverse, mixed_index_loss, Active, Active(index), Const(p), Const(x)
        )[1][1]
        step = cbrt(eps(T))
        expected =
            (mixed_index_loss(index + step, p, x) - mixed_index_loss(index - step, p, x)) /
            (2step)
        @test isfinite(actual)
        @test actual ≈ expected rtol=(T == Float32 ? 0.003 : 1e-7)
    end
end
