export PSNoiseGenerator

"""
    PSNoiseGenerator1D <: PowerSpectrumNoiseGenerator

This is a data type for a 1 dimensional power spectrum noise generator
"""

struct PSNoiseGenerator{P<:AbstractNoisePowerSpectrum, N<:NoiseSignal} <: PowerSpectrumNoiseGenerator
    psmodel::P
    noisesignal::N
    data
end

