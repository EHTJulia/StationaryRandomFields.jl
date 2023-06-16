export PSNoiseGenerator2D

"""
    PSNoiseGenerator2D <: PowerSpectrumNoiseGenerator

This is a data type for a 2 dimensional power spectrum noise generator
"""

struct PSNoiseGenerator2D <: PowerSpectrumNoiseGenerator
    psmodel::NoisePowerSpectrum
    noisesignal::NoiseSignal
    rfftplan 
    irfftplan
end

