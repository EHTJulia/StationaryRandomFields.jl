export AbstractContinuousNoiseSignal

"""
    AbstractContinuousNoiseSignal

This is an abstract data type for N dimensional continuous signal data on which signal noise will be generated. 
The data type extends the AbstractNoiseSignal data type for use in PowerSpectrumNoiseGenerator methods.

**Mandatory fields**
- `noisesignal::AbstractNoiseSignal`: noise signal object
"""

abstract type AbstractContinuousNoiseSignal end

###
### Methods to give the size
###
@inline function Base.size(data::AbstractContinuousNoiseSignal)::Tuple
    return data.dims
end

###
### Methods to give the size
###
@inline function Base.sizeof(data::AbstractContinuousNoiseSignal)::Integer
    return prod(data.dims)
end

###
### Methods to give the size
###
@inline function Base.ndims(data::AbstractContinuousNoiseSignal)::Integer
    return length(data.dims)
end

###
### Size of noises in Fourier domain
###
@inline function rfftsize(data::AbstractContinuousNoiseSignal)::Tuple
    return rfftsize(data.dims...)
end

###
### Compute frequencies in the Fourier Domain
###
@inline function AbstractFFTs.rfftfreq(data::AbstractContinuousNoiseSignal)::Tuple
    return rfftfreq(data.noisesignal)
end


###
### Generate a Complex Gaussian Noise
###
@inline function generate_gaussian_noise(data::AbstractContinuousNoiseSignal)
    return generate_gaussian_noise(data.noisesignal)
end 