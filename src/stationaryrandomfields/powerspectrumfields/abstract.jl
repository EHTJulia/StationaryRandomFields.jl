"""
    $(TYPEDEF)

**Mandatory fields**
- `dims::Tuple`: dimension of fields
- `dist::Distributions.UnivariateDistribution`: distribution of the random fields
- `fftplan`::plan for rfft
- `ifftplan`::inverse plan for rfft

This type is a random uncorrelated field where each element is drawn from a univariate distribution.
"""
abstract type AbstractPowerSpectrumRandomField{N} <: AbstractStationaryRandomField{N} end

@inline function Distributions._rand!(
    rng::AbstractRNG, field::AbstractPowerSpectrumRandomField, x::AbstractArray{<:Real}
)
    x[:] = forward(field, rand(rng, field.dist, field.dist.dims...))
    return nothing
end

@inline function Distributions._logpdf(
    field::AbstractPowerSpectrumRandomField, x::AbstractMatrix{<:Real}
)
    return Distributions._logpdf(field.dist, inverse(field, x))
end

@inline function forward(
    field::AbstractPowerSpectrumRandomField, xfourier::AbstractMatrix{<:Real}
)
    return field.ifftplan * xfourier
end

@inline function inverse(field::AbstractPowerSpectrumRandomField, x::AbstractMatrix{<:Real})
    return field.fftplan * x
end
