export NoiseSignal2D
export freq_norm

"""
    NoiseSignal2D <: AbstractNoiseSignal

This is a data type for two dimensional noise signals.


"""
struct NoiseSignal2D <: AbstractNoiseSignal
    dims::Tuple
    pixelsizes::Tuple
end

#
# Compute 2D grid of norms of each frequency point in the fourier plane
# 
@inline function freq_norm(data::NoiseSignal2D)::Tuple
    return ([u^2. + v^2. for u in rfftfreq(data)[1], v in rfftfreq(data)[2]],)
end

@inline function freq_norm(gridofν::Tuple)::Tuple
    return ([u^2. + v^2. for u in gridofν[1], v in gridofν[2]],)
end