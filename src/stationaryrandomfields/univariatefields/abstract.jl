export AbstractUnivariateRandomField

"""
    $(TYPEDEF)

This is an abstract field type for the `N`-dimensional stationary uncorrelated random fields where each element is drawn from a univariate distribution.
The field type contains the information for the dimension of the field and the distribution of the random fields.

**Mandatory fields**
- `dims::Tuple`: dimension of fields
- `dist::Distributions.UnivariateDistribution`: distribution of the random fields

**Mandatory methods**
See the documentation of `AbstractUnivariateRandomField`. No additional methods are required.
"""
abstract type AbstractUnivariateRandomField{N} <: AbstractStationaryRandomField{N} end

@inline function Distributions._rand!(
    rng::AbstractRNG, field::AbstractUnivariateRandomField, x::AbstractMatrix{<:Real}
)
    return Distributions.rand!(rng, field.dist, x)
end

@inline Distributions._logpdf(
    field::AbstractUnivariateRandomField, x::AbstractMatrix{<:Real}
) = sum(logpdf.(field.dist, x))