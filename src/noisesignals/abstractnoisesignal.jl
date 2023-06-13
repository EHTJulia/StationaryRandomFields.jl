export AbstractNoiseSignal
export generate_gaussian_noise
export rfftsize

"""
    AbstractNoiseSignal

This is an abstract data type for signal data on which powerspectrum noises will be created.
The data type contains information to define the dimensions and sizes of data, and provides
various methods.

# Mandatory fields
- `dims::Tuple`: dimension of data
- `pixelsizes::Tuple`: pixel size along with each dimension

# Mandatory methods
- `size(data::AbstractNoiseSignal)`: returns the size of the signal
- `ndims(data::AbstractNoiseSignal)`: returns the number of the dimension
- `xgrid(data::NoiseSignal1D)`: returns the position grid in the image plane
- `rfftsize(data::AbstractNoiseSignal)`: returns the size of the signal in the Fourier plane (for rfft)
- `rfftfreq(data::AbstractNoiseSignal)`: returns the freuqency grid along each dimension in the Fourier plane (for rfft)
- `generate_gaussian_noise(data::AbstractNoiseSignal)`: returns Gaussian noises in Fourier domain with the size of `rfftsize(data)`
"""
abstract type AbstractNoiseSignal end

###
### Methods to give the size
###
@inline function Base.size(data::AbstractNoiseSignal)::Tuple
    return data.dims
end

###
### Methods to give the size
###
#@inline function Base.sizeof(data::AbstractNoiseSignal)::Tuple
#    return prod(data.dims)
#end

###
### Methods to give the size
###
@inline function Base.ndims(data::AbstractNoiseSignal)::Integer
    return length(data.dims)
end

#
# Generate x-axis grid in position space for graphing data 
#
@inline function xgrid(data::AbstractNoiseSignal)
    return [i*data.pixelsizes for i in 0:(data.dims-1)]
end 

###
### Size of noises in Fourier domain
###
"""
    rfftsize(data::AbstractNoiseSignal)::Tuple

Returns the size of the Fourier-domain signal
"""
@inline function rfftsize(data::AbstractNoiseSignal)::Tuple
    return rfftsize(data.dims...)
end

@inline function rfftsize(dims::Number)::Tuple
    return tuple(dims ÷ 2 + 1,)
end

@inline function rfftsize(dims::Number...)::Tuple
    if length(dims) > 1
        return (dims[1] ÷ 2 + 1, dims[2:end]...)
    end
end

###
### Compute frequencies in the Fourier Domain
###
@inline function AbstractFFTs.rfftfreq(data::AbstractNoiseSignal)::Tuple
    return Tuple((
        if i == 1
            rfftfreq(data.dims[i], data.pixelsizes[i])
        else
            fftfreq(data.dims[i], data.pixelsizes[i])
        end
    ) for i in 1:ndims(data))
end

###
### Generate a Complex Gaussian Noise
###
"""
    generate_gaussian_noise(data::AbstractNoiseSignal; rng=Random.default_rng())

Generate complex gausisian noises on Fourier space over the frequency space expected for Real FFT.
This return a complex array, which is consistent with Fourier-domain representation of a real uncorrelated Gaussian noise
with the zero mean and variance of `sizeof(data)/2` in signal-domain. The size of the output array is given by `rfftsize(data)`.
"""
@inline function generate_gaussian_noise(data::AbstractNoiseSignal; rng=Random.default_rng())
    # generate Gaussian noises
    ϵre = randn(rng, Float64, rfftsize(data)...)
    ϵim = randn(rng, Float64, rfftsize(data)...)
    ϵre[1] = 0.0 # Mean of the generated signals will be 0
    ϵim[1] = 0.0 # Real signal should be real at this frequency

    # Fill zero at frequencies where fourier component must be real for real signals
    evenness = iseven.(data.dims)
    midindex = Tuple(dim ÷ 2 + 1 for dim in data.dims)

    if evenness[1]
        # create iterator
        iter = []
        for idim in 1:ndims(data)
            if evenness[idim]
                push!(iter, tuple(1, midindex[idim]))
            else
                push!(iter, tuple(1,))
            end
        end

        # fill zero
        for idx in Iterators.product(iter...)
            ϵim[idx...] = 0.0
        end
    end

    return ϵre .+ im .* ϵim
end
