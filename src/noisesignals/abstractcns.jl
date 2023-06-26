export AbstractContinuousNoiseSignal

"""
    AbstractContinuousNoiseSignal

This is an abstract data type for N dimensional continuous signal data on which signal noise will be generated. 
The data type extends the AbstractNoiseSignal data type for use in PowerSpectrumNoiseGenerator methods.

**Mandatory fields**
- `noisesignal::AbstractNoiseSignal`: noise signal object
"""

abstract type AbstractContinuousNoiseSignal end


@inline function Base.size(data::AbstractContinuousNoiseSignal)::Tuple
    return data.dims
end


@inline function Base.sizeof(data::AbstractContinuousNoiseSignal)::Integer
    return prod(data.dims)
end


@inline function Base.ndims(data::AbstractContinuousNoiseSignal)::Integer
    return length(data.dims)
end


@inline function rfftsize(data::AbstractContinuousNoiseSignal)::Tuple
    return rfftsize(data.dims...)
end


@inline function AbstractFFTs.rfftfreq(data::AbstractContinuousNoiseSignal)::Tuple
    return rfftfreq(data.noisesignal)
end


@inline function generate_gaussian_noise(data::AbstractContinuousNoiseSignal)
    return generate_gaussian_noise(data.noisesignal)
end 