export NoiseSignal1D

"""
    NoiseSignal1D <: AbstractNoiseSignal

This is a data type for N dimensional noise signals.
"""
struct NoiseSignal <: AbstractNoiseSignal
    dims::Tuple
end