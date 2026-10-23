"""An array-valued continuous distribution with one sample of rank `N`."""
abstract type AbstractStationaryRandomField{N} <:
              Distribution{ArrayLikeVariate{N},Continuous} end

"""
    UnivariateRandomField(dims, dist)

An independent array of identically distributed continuous scalar variates.
`Univariate` describes each element, not the joint distribution of a sample.
"""
struct UnivariateRandomField{N,D<:ContinuousUnivariateDistribution} <:
       AbstractStationaryRandomField{N}
    dims::NTuple{N,Int}
    dist::D
    function UnivariateRandomField(
        dims::NTuple{N,Integer}, dist::D
    ) where {N,D<:ContinuousUnivariateDistribution}
        N > 0 || throw(ArgumentError("A field must have at least one dimension."))
        all(>(0), dims) || throw(ArgumentError("Field dimensions must be positive."))
        return new{N,D}(Int.(dims), dist)
    end
end

Base.size(f::UnivariateRandomField) = f.dims
function Base.size(f::AbstractStationaryRandomField, d::Integer)
    d > 0 ? (d <= ndims(f) ? size(f)[d] : 1) : throw(ArgumentError("Invalid axis."))
end
Base.ndims(::AbstractStationaryRandomField{N}) where {N} = N
Base.length(f::AbstractStationaryRandomField) = prod(size(f))
distribution_eltype(::Type{D}) where {D<:ContinuousUnivariateDistribution} = eltype(D)
distribution_eltype(::Type{Normal{T}}) where {T} = T
distribution_eltype(::Type{Laplace{T}}) where {T} = T
Base.eltype(::Type{UnivariateRandomField{N,D}}) where {N,D} = distribution_eltype(D)
Base.eltype(f::UnivariateRandomField) = eltype(typeof(f))

"""Construct an independent Gaussian field. `σ` is a standard deviation."""
function UnivariateGaussianRandomField(dims; μ=0.0, σ=1.0)
    UnivariateRandomField(dims, Normal(μ, σ))
end
function UnivariateGaussianRandomField(
    ::Type{T}, dims; μ=zero(T), σ=one(T)
) where {T<:AbstractFloat}
    UnivariateRandomField(dims, Normal(T(μ), T(σ)))
end
UnivariateGaussianRandomField(dims, d::Normal) = UnivariateRandomField(dims, d)

"""Construct an independent Laplace field. `b` is a scale, with variance `2b^2`."""
function UnivariateLaplaceRandomField(dims; μ=0.0, b=1.0)
    UnivariateRandomField(dims, Laplace(μ, b))
end
function UnivariateLaplaceRandomField(
    ::Type{T}, dims; μ=zero(T), b=one(T)
) where {T<:AbstractFloat}
    UnivariateRandomField(dims, Laplace(T(μ), T(b)))
end
UnivariateLaplaceRandomField(dims, d::Laplace) = UnivariateRandomField(dims, d)

function Distributions._rand!(
    rng::AbstractRNG, f::UnivariateRandomField{N}, x::AbstractArray{<:Real,N}
) where {N}
    rand!(rng, f.dist, x)
    return x
end

# The scalar formulas remain visible to device compilers and differentiation.
function latent_logpdf(d::Normal, x)
    μ, σ = mean(d), std(d)
    return -sum(v -> abs2((v - μ) / σ), x) / 2 - length(x) * log(σ * sqrt(oftype(σ, 2π)))
end
latent_logpdf(d::Laplace, x) = -sum(v -> abs(v - d.μ), x) / d.θ - length(x) * log(2 * d.θ)
latent_logpdf(d::ContinuousUnivariateDistribution, x) = sum(v -> logpdf(d, v), x)

function Distributions._logpdf(
    f::UnivariateRandomField{N}, x::AbstractArray{<:Real,N}
) where {N}
    latent_logpdf(f.dist, x)
end
# Distributions has a matrix-specific fallback, requiring this disambiguation.
function Distributions._logpdf(f::UnivariateRandomField{2}, x::AbstractMatrix{<:Real})
    latent_logpdf(f.dist, x)
end

function Distributions.insupport(
    f::UnivariateRandomField{N}, x::AbstractArray{<:Real,N}
) where {N}
    return size(x) == size(f) && all(v -> insupport(f.dist, v), x)
end
Distributions.insupport(::UnivariateRandomField, ::Any) = false

Base.zeros(f::AbstractStationaryRandomField) = zeros(eltype(f), size(f))
Base.ones(f::AbstractStationaryRandomField) = ones(eltype(f), size(f))
Base.fill(v, f::AbstractStationaryRandomField) = fill(convert(eltype(f), v), size(f))
