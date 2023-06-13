export PSNoiseGenerator1D
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

@inline function generate_signal_noise(psgen::PSNoiseGenerator1D)
    return psgen.irfftplan * (map_ampspectrum(psgen.psmodel, psgen.noisesignal)[1] .* generate_gaussian_noise(psgen.noisesignal))
end