export PSNoiseGenerator1D


"""
    PSNoiseGenerator1D <: PowerSpectrumNoiseGenerator

This is a data type for a 1 dimensional power spectrum noise generator
"""

struct PSNoiseGenerator1D <: PowerSpectrumNoiseGenerator
    psmodel::AbstractNoisePowerSpectrum
    noisesignal::NoiseSignal
    rfftplan 
    irfftplan
end

