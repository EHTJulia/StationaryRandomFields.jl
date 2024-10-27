export AbstractStationaryRandomField
export length, eltype, _rand!, _logpdf
export size, sizeof, ndims
export zeros, ones, Array, fill

"""
    AbstractStationaryRandomField{N} <: Distributions.Distribution{Distributions.ArrayLikeVariate{N}, Distributions.Continuous}

This is an abstract field type for the `N`-dimensional stationary random fields.
The field type contains information to define the size and dimension of the field.

**Mandatory fields**
- `dims::Tuple`: dimension of fields

**Mandatory methods**
As a subtype of `Distributions.Distribution`, it should have the following methods.
- `length(d::AbstractStationaryRandomField)`
- `eltype(d::AbstractStationaryRandomField)`
- `Distributions._rand!(::AbstractRNG, d::AbstractStationaryRandomField, x::AbstractArray)``
- `Distributions._logpdf(d::AbstractStationaryRandomField, x::AbstractArray)``

Here are the additional methods that should be implemented.
- `size`: returns the size of the fields
- `sizeof`: returns product of signal dimensions
- `ndims`: returns the number of the dimension
"""
abstract type AbstractStationaryRandomField{N} <: Distributions.Distribution{Distributions.ArrayLikeVariate{N}, Distributions.Continuous} end

#
# Mandatory methods for Distributions.Distribution
#
@inline Base.length(::AbstractStationaryRandomField{N}) where N = N

# Each subtype should implement this method
function Base.eltype(::AbstractStationaryRandomField) end

# Each subtype should implement this method
#function Distributions._rand!(::AbstractRNG, ::AbstractStationaryRandomField, ::AbstractArray) end

# Each subtype should implement this method
#function Distributions._logpdf!(::AbstractRNG, ::AbstractStationaryRandomField, ::AbstractArray) end

#
# Mandatory methods for AbstractStationaryRandomField
#
@inline Base.size(field::AbstractStationaryRandomField) = field.dims

@inline Base.sizeof(field::AbstractStationaryRandomField) = prod(field.dims)

@inline Base.ndims(field::AbstractStationaryRandomField) = length(field)

#
# Other useful methods
#
@inline Base.zeros(field::AbstractStationaryRandomField) = zeros(eltype(field), field.dims...)
@inline Base.ones(field::AbstractStationaryRandomField) = ones(eltype(field), field.dims...)
@inline Base.Array(undef, field) = Array{eltype(field)}(undef, field.dims...)
@inline Base.fill(value, field::AbstractStationaryRandomField) = fill(eltype(field)(value), field.dims...)