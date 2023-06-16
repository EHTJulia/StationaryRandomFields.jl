export NoiseSignal2D
export freq_norm

"""
    NoiseSignal2D <: AbstractNoiseSignal

This is a data type for two dimensional noise signals.

# Mandatory methods
- `freq_norm`: 2D matrix of norms of each frequency point in fourier space

"""
struct NoiseSignal2D <: AbstractNoiseSignal
    dims::Tuple
    pixelsizes::Tuple
end

#
# Compute grid of norms of each frequency point in the fourier plane (currently only works for 1D and 2D)
# 
@inline function freq_norm(data::NoiseSignal2D)::Array
    return freq_norm(rfftfreq(data))
end

@inline function freq_norm(gridofν::Tuple)::Array
    if length(gridofν)==1
        return gridofν
    else
        return [u^2. + v^2. for u in gridofν[1], v in gridofν[2]]
    end
end

