export AbstractContinuousNoiseSignal

"""
AbstractContinuousNoiseSignal

This is an abstract data type for N dimensional continuous signal data on which signal noise will be generated. 
The data type contains an AbstractNoiseSignal and defines its corresponding dimensions, RFFT plan, and inverse RFFT plan. 
The type of data (default Float64) may also be specified.

# Mandatory fields
- `noisesignal::AbstractNoiseSignal`: previously defined noise signal object
- `dims::Tuple`: dimensions of noisesignal
- `plan`: rfft plan corresponding to noisesignal
- `invplan`: irfft plan corresponding to noisesignal 

# Mandatory methods
- `size(data::AbstractContinuousNoiseSignal)`: returns the size of the signal
- `sizeof(data::AbstractContinuousNoiseSignal)`: returns product of signal dimensions
- `ndims(data::AbstractContinuousNoiseSignal)`: returns the number of the dimension
- `rfftsize(data::AbstractContinuousNoiseSignal)`: returns the size of the signal in the Fourier plane (for rfft)
- `rfftfreq(data::AbstractContinuousNoiseSignal)`: returns the frequency grid along each dimension in the Fourier plane (for rfft)
- `generate_gaussian_noise(data::AbstractContinuousNoiseSignal)`: returns Gaussian noises in Fourier domain with the size of `rfftsize(data)`

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

#
# Compute grid of norms of each frequency point in the fourier plane (currently only works for 1D and 2D)
# 
@inline function freq_norm(data::AbstractContinuousNoiseSignal)::Array
    return freq_norm(data.noisesignal)
end

###
### Generate a Complex Gaussian Noise
###

@inline function generate_gaussian_noise(data::AbstractContinuousNoiseSignal)
    return generate_gaussian_noise(data.noisesignal)
end 