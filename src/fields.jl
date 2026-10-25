"""An array-valued continuous distribution with one sample of rank `N`."""
abstract type AbstractStationaryRandomField{N} <:
Distribution{ArrayLikeVariate{N}, Continuous} end

"""
    UnivariateRandomField(dims, dist)

An independent array of identically distributed continuous scalar variates.
`Univariate` describes each element, not the joint distribution of a sample.
"""
struct UnivariateRandomField{N, D <: ContinuousUnivariateDistribution} <:
    AbstractStationaryRandomField{N}
    dims::NTuple{N, Int}
    dist::D
    function UnivariateRandomField(
            dims::NTuple{N, Integer}, dist::D
        ) where {N, D <: ContinuousUnivariateDistribution}
        N > 0 || throw(ArgumentError("A field must have at least one dimension."))
        all(>(0), dims) || throw(ArgumentError("Field dimensions must be positive."))
        return new{N, D}(Int.(dims), dist)
    end
end

Base.size(f::UnivariateRandomField) = f.dims
function Base.size(f::AbstractStationaryRandomField, d::Integer)
    return d > 0 ? (d <= ndims(f) ? size(f)[d] : 1) : throw(ArgumentError("Invalid axis."))
end
Base.ndims(::AbstractStationaryRandomField{N}) where {N} = N
Base.length(f::AbstractStationaryRandomField) = prod(size(f))
distribution_eltype(::Type{D}) where {D <: ContinuousUnivariateDistribution} = eltype(D)
distribution_eltype(::Type{Normal{T}}) where {T} = T
distribution_eltype(::Type{Laplace{T}}) where {T} = T
Base.eltype(::Type{UnivariateRandomField{N, D}}) where {N, D} = distribution_eltype(D)
Base.eltype(f::UnivariateRandomField) = eltype(typeof(f))

"""Construct an independent Gaussian field. `σ` is a standard deviation."""
function UnivariateGaussianRandomField(dims; μ = 0.0, σ = 1.0)
    return UnivariateRandomField(dims, Normal(μ, σ))
end
function UnivariateGaussianRandomField(
        ::Type{T}, dims; μ = zero(T), σ = one(T)
    ) where {T <: AbstractFloat}
    return UnivariateRandomField(dims, Normal(T(μ), T(σ)))
end
UnivariateGaussianRandomField(dims, d::Normal) = UnivariateRandomField(dims, d)

"""Construct an independent Laplace field. `b` is a scale, with variance `2b^2`."""
function UnivariateLaplaceRandomField(dims; μ = 0.0, b = 1.0)
    return UnivariateRandomField(dims, Laplace(μ, b))
end
function UnivariateLaplaceRandomField(
        ::Type{T}, dims; μ = zero(T), b = one(T)
    ) where {T <: AbstractFloat}
    return UnivariateRandomField(dims, Laplace(T(μ), T(b)))
end
UnivariateLaplaceRandomField(dims, d::Laplace) = UnivariateRandomField(dims, d)

function Distributions._rand!(
        rng::AbstractRNG, f::UnivariateRandomField{N}, x::AbstractArray{<:Real, N}
    ) where {N}
    rand!(rng, f.dist, x)
    return x
end
function _rand_normal_field!(rng, f, x)
    Random.randn!(rng, x)
    x .= mean(f.dist) .+ std(f.dist) .* x
    return x
end
function _rand_laplace_field!(rng, f, x)
    Random.rand!(rng, x)
    # Map the midpoint of each representable uniform bin into (0, 1).
    x .= clamp.(x, eps(eltype(x)), one(eltype(x)) - eps(eltype(x))) .- 0.5
    x .= f.dist.μ .- f.dist.θ .* sign.(x) .* log1p.(-2 .* abs.(x))
    return x
end

# The scalar formulas remain visible to device compilers and differentiation.
function latent_logpdf(d::Normal, x)
    return latent_logpdf(d, x, mean(d), std(d))
end
function latent_logpdf(::Normal, x, μ, σ)
    return -sum(v -> abs2((v - μ) / σ), x) / 2 - length(x) * log(σ * sqrt(oftype(σ, 2π)))
end
latent_logpdf(d::Laplace, x) = latent_logpdf(d, x, d.μ, d.θ)
latent_logpdf(::Laplace, x, μ, b) = -sum(v -> abs(v - μ), x) / b - length(x) * log(2 * b)
prior_scale(d::Normal) = std(d)
prior_scale(d::Laplace) = d.θ

function parameterized_logpdf(f, x, location, scale, check_args)
    size(x) == size(f) || throw(DimensionMismatch("Wrong latent shape."))
    (!check_args || (isfinite(scale) && scale > 0 && isfinite(location))) ||
        throw(ArgumentError("Require a finite location and a positive finite scale."))
    return latent_logpdf(f.dist, x, location, scale)
end

"""
    logpdf(field, x; location = mean(field.dist), scale, check_args = true)

Evaluate an independent Gaussian or Laplace latent prior with optional active
location and scale parameters. `scale` is Gaussian standard deviation or Laplace
scale and defaults to the stored value. This avoids constructing a Distributions
object from traced scalar parameters. Set `check_args=false` inside a compiled
objective after enforcing a positive scale, for example with `exp(logscale)`.
"""
function Distributions.logpdf(
        f::UnivariateRandomField{N, D},
        x::AbstractArray{<:Real, N};
        location = mean(f.dist),
        scale = prior_scale(f.dist),
        check_args = true,
    ) where {N, D <: Union{Normal, Laplace}}
    return parameterized_logpdf(f, x, location, scale, check_args)
end
latent_logpdf(d::ContinuousUnivariateDistribution, x) = sum(v -> logpdf(d, v), x)

function Distributions._logpdf(
        f::UnivariateRandomField{N}, x::AbstractArray{<:Real, N}
    ) where {N}
    return latent_logpdf(f.dist, x)
end
# Distributions has a matrix-specific fallback, requiring this disambiguation.
function Distributions._logpdf(f::UnivariateRandomField{2}, x::AbstractMatrix{<:Real})
    return latent_logpdf(f.dist, x)
end

function Distributions.insupport(
        f::UnivariateRandomField{N}, x::AbstractArray{<:Real, N}
    ) where {N}
    return size(x) == size(f) && all(v -> insupport(f.dist, v), x)
end
Distributions.insupport(::UnivariateRandomField, ::Any) = false

Base.zeros(f::AbstractStationaryRandomField) = zeros(eltype(f), size(f))
Base.ones(f::AbstractStationaryRandomField) = ones(eltype(f), size(f))
Base.fill(v, f::AbstractStationaryRandomField) = fill(convert(eltype(f), v), size(f))
