using Enzyme

function field_loss(x, h, s, p, target)
    return sum(abs2, apply(h, apply(s, x, p)) .- target)
end
function finite_gradient(f, x)
    h = cbrt(eps(eltype(x)))
    return (v -> reshape(v, size(x)))(
        map(eachindex(x)) do i
            a, b = copy(x), copy(x)
            a[i] += h
            b[i] -= h
            (f(a) - f(b)) / (2h)
        end,
    )
end
function latent_loss(x, f)
    return logpdf(f, x)
end
function wavelet_loss(x, p)
    return sum(abs2, apply(p, x) .- 0.3)
end
function spectrum_loss(a, x, h, s)
    model = ConstantPowerSpectrum{ndims(x)}(exp(a))
    return sum(abs2, apply(h, apply(s, x, model)) .- 0.3)
end

@testset "Enzyme latent and transform gradients" begin
    for T in (Float32, Float64)
        x = reshape(T.(0.2 .+ sin.(1:16)), 4, 4)
        for dist in (Normal(T(0.1), T(1.3)), Laplace(T(-0.2), T(0.7)))
            f = UnivariateRandomField(size(x), dist)
            dx = zero(x)
            Enzyme.autodiff(Reverse, latent_loss, Active, Duplicated(x, dx), Const(f))
            expected = if dist isa Normal
                -(x .- mean(dist)) ./ var(dist)
            else
                -sign.(x .- dist.μ) ./ dist.θ
            end
            @test dx ≈ expected
        end
        h = prepare(HartleyTransform(size(x)), x)
        s = prepare(SpectralScaling(FieldGrid(size(x)), ConstantPowerSpectrum{2}(T(2))), x)
        target = fill(T(0.1), size(x))
        dx = zero(x)
        p = ConstantPowerSpectrum{2}(T(2))
        Enzyme.autodiff(
            Reverse,
            field_loss,
            Active,
            Duplicated(x, dx),
            Const(h),
            Const(s),
            Const(p),
            Const(target),
        )
        @test dx ≈ 2sqrt(T(2)) .* apply(h, sqrt(T(2)) .* apply(h, x) .- target)
        tangent = Enzyme.autodiff(
            Forward,
            field_loss,
            Duplicated,
            Duplicated(x, one.(x)),
            Const(h),
            Const(s),
            Const(p),
            Const(target),
        )[1]
        @test tangent ≈ sum(dx)
        for wt in (wavelet(WT.haar), wavelet(WT.db2))
            w = prepare(WaveletTransform(size(x), wt; levels=2), x)
            dx = zero(x)
            Enzyme.autodiff(Reverse, wavelet_loss, Active, Duplicated(x, dx), Const(w))
            @test dx ≈ SRF.wavelet_analysis(w, 2 .* (apply(w, x) .- 0.3))
            tangent = Enzyme.autodiff(
                Forward, wavelet_loss, Duplicated, Duplicated(x, one.(x)), Const(w)
            )[1]
            @test tangent ≈ sum(dx)
        end
        a = T(0.4)
        da = Enzyme.autodiff(
            Reverse, spectrum_loss, Active, Active(a), Const(x), Const(h), Const(s)
        )[1][1]
        step = cbrt(eps(T))
        finite =
            (spectrum_loss(a + step, x, h, s) - spectrum_loss(a - step, x, h, s)) / (2step)
        @test isapprox(da, finite; rtol=T == Float32 ? 0.003 : 1e-7)
    end
end

function indexed_spectrum_loss(a, x, h, s)
    model = SinglePowerLaw{ndims(x)}(a; dc=0.0)
    return sum(abs2, apply(h, apply(s, x, model)) .- 0.3)
end
@testset "Spectrum index gradient with a zero DC mode" begin
    x = reshape(sin.(1:20), 4, 5)
    h = prepare(HartleyTransform(size(x)), x)
    s = prepare(SpectralScaling(FieldGrid(size(x)), SinglePowerLaw{2}(-1.2)), x)
    a = -1.2
    da = Enzyme.autodiff(
        Reverse, indexed_spectrum_loss, Active, Active(a), Const(x), Const(h), Const(s)
    )[1][1]
    step = 1e-5
    finite =
        (
            indexed_spectrum_loss(a + step, x, h, s) -
            indexed_spectrum_loss(a - step, x, h, s)
        ) / (2step)
    @test isfinite(da)
    @test da ≈ finite rtol=1e-7
end

function geometry_loss(a, x, h, s)
    model = modify(
        SinglePowerLaw{2}(-1.2),
        Stretch(exp(sum(a[1:1])), exp(sum(a[2:2])); check_args=false),
        Rotate(sum(a[3:3])),
    )
    return sum(abs2, apply(h, apply(s, x, model)) .- 0.3)
end
@testset "Spectrum stretch and angle gradients" begin
    x = reshape(sin.(1:20), 4, 5)
    h = prepare(HartleyTransform(size(x)), x)
    s = prepare(SpectralScaling(FieldGrid(size(x)), SinglePowerLaw{2}(-1.2)), x)
    a = [0.1, -0.2, 0.3]
    da = zero(a)
    Enzyme.autodiff(
        Reverse, geometry_loss, Active, Duplicated(a, da), Const(x), Const(h), Const(s)
    )
    expected = finite_gradient(b -> geometry_loss(b, x, h, s), a)
    @test all(isfinite, da)
    @test da ≈ expected rtol=1e-7
end

function prior_parameter_loss(a, x, ::Val{family}) where {family}
    μ, scale = sum(a[1:1]), exp(sum(a[2:2]))
    dist = family === :normal ? Normal() : Laplace()
    field = UnivariateRandomField(size(x), dist)
    return logpdf(field, x; location=μ, scale=scale, check_args=false)
end
function prior_parameter_gradient(a, x, family)
    da = zero(a)
    Enzyme.autodiff(
        Reverse, prior_parameter_loss, Active, Duplicated(a, da), Const(x), Const(family)
    )
    return da
end
@testset "Latent prior location and log-scale gradients" begin
    for T in (Float32, Float64), family in (Val(:normal), Val(:laplace))
        x = reshape(T.(0.2 .+ sin.(1:16)), 4, 4)
        a = T[0.1, -0.2]
        actual = prior_parameter_gradient(a, x, family)
        expected = finite_gradient(b -> prior_parameter_loss(b, x, family), a)
        @test all(isfinite, actual)
        @test actual ≈ expected rtol = (T == Float32 ? 0.003 : 1e-7)
    end
end
