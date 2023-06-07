export NoiseSignal2D

"""
    NoiseSignal2D <: AbstractNoiseSignal

This is a data type for one dimensional signals 
"""
struct NoiseSignal2D <: AbstractNoiseSignal
    dims::Tuple
    pixelsizes::Tuple
end