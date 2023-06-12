"""
    PSNoiseGenerator1D <: PowerSpectrumNoiseGenerator

This is a data type for a 1 dimensional power spectrum noise generator
"""

struct PSNoiseGenerator1D <: PowerSpectrumNoiseGenerator
    psmodel::AbstractNoisePowerSpectrum
    noisesignal::AbstractNoiseSignal
    rfftplan 
    irfftplan
end