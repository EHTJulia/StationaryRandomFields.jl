export NoiseSignal

"""
    NoiseSignal <: AbstractNoiseSignal

This is a data type for N dimensional noise signals.
"""
struct NoiseSignal <: AbstractNoiseSignal
    dims::Tuple
end