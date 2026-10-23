using Test, Reactant, Enzyme, Distributions, Random, Wavelets, StationaryRandomFields

Reactant.set_default_backend(get(ENV, "SRF_REACTANT_BACKEND", "cpu"))
const SRF = StationaryRandomFields

function latent_gradient(x, f)
    dx = zero(x)
    Enzyme.autodiff(Reverse, logpdf, Active, Const(f), Duplicated(x, dx))
    return dx
end
function transform_loss(x, p)
    return sum(abs2, apply(p, x) .- 0.3)
end
function transform_gradient(x, p)
    dx = zero(x)
    Enzyme.autodiff(Reverse, transform_loss, Active, Duplicated(x, dx), Const(p))
    return dx
end
function parameter_loss(a, x, h, s)
    return sum(
        abs2,
        apply(
            h, apply(s, x, ConstantPowerSpectrum{ndims(x)}(exp(sum(a)); check_args=false))
        ) .- 0.3,
    )
end
function parameter_gradient(a, x, h, s)
    da = zero(a)
    Enzyme.autodiff(
        Reverse, parameter_loss, Active, Duplicated(a, da), Const(x), Const(h), Const(s)
    )
    return da
end
function draw!(rng, f, x)
    return rand!(rng, f, x)
end

@testset "Reactant CPU/GPU field operations" begin
    println(
        "Julia=",
        VERSION,
        " Reactant=",
        pkgversion(Reactant),
        " Enzyme=",
        pkgversion(Enzyme),
        " backend=",
        get(ENV, "SRF_REACTANT_BACKEND", "cpu"),
    )
    for T in (Float32, Float64), dims in ((8,), (4, 4), (4, 2, 2), (2, 2, 2, 2))
        x = reshape(T.(0.2 .+ sin.(1:prod(dims))), dims)
        rx = Reactant.to_rarray(x)
        for dist in (Normal(T(0.1), T(1.3)), Laplace(T(-0.2), T(0.7)))
            f = UnivariateRandomField(dims, dist)
            compiled = Reactant.@compile logpdf(f, rx)
            @test T(compiled(f, rx)) ≈ logpdf(f, x)
            gradient = Reactant.@compile latent_gradient(rx, f)
            expected = if dist isa Normal
                -(x .- mean(dist)) ./ var(dist)
            else
                -sign.(x .- dist.μ) ./ dist.θ
            end
            @test Array(gradient(rx, f)) ≈ expected
        end
        for t in (HartleyTransform(dims), WaveletTransform(dims; levels=1))
            p = prepare(t, rx)
            compiled = Reactant.@compile apply(p, rx)
            expected = apply(prepare(t, x), x)
            @test Array(compiled(p, rx)) ≈ expected
            gradient = Reactant.@compile transform_gradient(rx, p)
            analytical = if t isa HartleyTransform
                apply(prepare(t, x), 2 .* (expected .- T(0.3)))
            else
                SRF.wavelet_analysis(prepare(t, x), 2 .* (expected .- T(0.3)))
            end
            @test Array(gradient(rx, p)) ≈ analytical
        end
        h = prepare(HartleyTransform(dims), rx)
        s = prepare(
            SpectralScaling(FieldGrid(dims), ConstantPowerSpectrum{length(dims)}(T(2))), rx
        )
        a = Reactant.to_rarray(T[0.4])
        gradient = Reactant.@compile parameter_gradient(a, rx, h, s)
        result = only(Array(gradient(a, rx, h, s)))
        hc = prepare(HartleyTransform(dims), x)
        sc = prepare(
            SpectralScaling(FieldGrid(dims), ConstantPowerSpectrum{length(dims)}(T(2))), x
        )
        step = cbrt(eps(T))
        finite =
            (
                parameter_loss(T[0.4 + step], x, hc, sc) -
                parameter_loss(T[0.4 - step], x, hc, sc)
            ) / (2step)
        @test isapprox(result, finite; rtol=T == Float32 ? 0.003 : 1e-7)
    end
    for dist in (Normal(0.0f0, 1.0f0), Laplace(0.0f0, 1.0f0))
        f = UnivariateRandomField((64,), dist)
        rng = Reactant.ReactantRNG(Reactant.to_rarray(UInt64[12, 34]))
        x = Reactant.to_rarray(zeros(Float32, 64))
        draw = Reactant.@compile draw!(rng, f, x)
        draw(rng, f, x)
        first = Array(x)
        draw(rng, f, x)
        @test all(isfinite, Array(x))
        @test first != Array(x)
    end
end

function index_loss(a, x, h, s)
    model = SinglePowerLaw{ndims(x)}(sum(a); dc=0.0, check_args=false)
    return sum(abs2, apply(h, apply(s, x, model)) .- 0.3)
end
function index_gradient(a, x, h, s)
    da = zero(a)
    Enzyme.autodiff(
        Reverse, index_loss, Active, Duplicated(a, da), Const(x), Const(h), Const(s)
    )
    return da
end
@testset "Compiled power-law index gradient" begin
    x = reshape(sin.(1:20), 4, 5)
    rx = Reactant.to_rarray(x)
    h = prepare(HartleyTransform(size(x)), rx)
    s = prepare(SpectralScaling(FieldGrid(size(x)), SinglePowerLaw{2}(-1.2)), rx)
    a = Reactant.to_rarray([-1.2])
    gradient = Reactant.@compile index_gradient(a, rx, h, s)
    result = only(Array(gradient(a, rx, h, s)))
    hc = prepare(HartleyTransform(size(x)), x)
    sc = prepare(SpectralScaling(FieldGrid(size(x)), SinglePowerLaw{2}(-1.2)), x)
    step = 1e-5
    expected =
        (index_loss([-1.2 + step], x, hc, sc) - index_loss([-1.2 - step], x, hc, sc)) /
        (2step)
    @test isfinite(result)
    @test result ≈ expected rtol=1e-7
end

function buffered(p, x)
    y = similar(x)
    return apply!(y, p, x)
end
@testset "Compiled buffered synthesis and signal sampling" begin
    x = reshape(Float32.(sin.(1:16)), 4, 4)
    rx = Reactant.to_rarray(x)
    for t in (HartleyTransform(size(x)), WaveletTransform(size(x); levels=2))
        p = prepare(t, rx)
        compiled = Reactant.@compile buffered(p, rx)
        @test Array(compiled(p, rx)) ≈ apply(prepare(t, x), x)
    end
    latent = UnivariateGaussianRandomField(Float32, size(x))
    descriptor = TransformChain(
        SpectralScaling(FieldGrid(size(x)), ConstantPowerSpectrum{2}(1.0f0)),
        HartleyTransform(size(x)),
    )
    field = prepare(TransformedRandomField(latent, descriptor), rx)
    rng = Reactant.ReactantRNG(Reactant.to_rarray(UInt64[12, 34]))
    draw = Reactant.@compile draw!(rng, field, rx)
    draw(rng, field, rx)
    first = Array(rx)
    draw(rng, field, rx)
    @test all(isfinite, Array(rx))
    @test first != Array(rx)
end

function geometry_loss(a, x, h, s)
    model = modify(
        SinglePowerLaw{2}(-1.2),
        Stretch(exp(sum(a[1:1])), exp(sum(a[2:2])); check_args=false),
        Rotate(sum(a[3:3])),
    )
    return sum(abs2, apply(h, apply(s, x, model)) .- 0.3)
end
function geometry_gradient(a, x, h, s)
    da = zero(a)
    Enzyme.autodiff(
        Reverse, geometry_loss, Active, Duplicated(a, da), Const(x), Const(h), Const(s)
    )
    return da
end
@testset "Compiled spectrum stretch and angle gradients" begin
    x = reshape(sin.(1:20), 4, 5)
    rx = Reactant.to_rarray(x)
    h = prepare(HartleyTransform(size(x)), rx)
    s = prepare(SpectralScaling(FieldGrid(size(x)), SinglePowerLaw{2}(-1.2)), rx)
    a = [0.1, -0.2, 0.3]
    ra = Reactant.to_rarray(a)
    gradient = Reactant.@compile geometry_gradient(ra, rx, h, s)
    actual = Array(gradient(ra, rx, h, s))
    hc = prepare(HartleyTransform(size(x)), x)
    sc = prepare(SpectralScaling(FieldGrid(size(x)), SinglePowerLaw{2}(-1.2)), x)
    step = 1e-5
    expected = map(eachindex(a)) do i
        plus, minus = copy(a), copy(a)
        plus[i] += step
        minus[i] -= step
        (geometry_loss(plus, x, hc, sc) - geometry_loss(minus, x, hc, sc)) / (2step)
    end
    @test all(isfinite, actual)
    @test actual ≈ expected rtol=1e-7
end
