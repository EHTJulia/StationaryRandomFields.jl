export PSNoiseGenerator

"""
    PSNoiseGenerator1D <: PowerSpectrumNoiseGenerator

This is a data type for a 1 dimensional power spectrum noise generator
"""

struct PSNoiseGenerator{P<:AbstractPowerSpectrumModel, N<:ContinuousNoiseSignal} <: PowerSpectrumNoiseGenerator
    psmodel::P
    noisesignal::N
end

