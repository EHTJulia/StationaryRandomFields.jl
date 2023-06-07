export NoiseSignal1D

"""
    NoiseSignal1D <: AbstractNoiseSignal

This is a data type for one dimensional signals 
"""
struct NoiseSignal1D <: AbstractNoiseSignal
    dims::Tuple
    pixelsizes::Tuple
end