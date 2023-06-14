export PSNoiseGenerator1D
export generate_signal_noise
export get_power_spectrum

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

#
# Convert fourier domain gaussian noise to signal noise with power law and inverse rfft 
#
@inline function generate_signal_noise(psgen::PSNoiseGenerator1D)::Vector
    return psgen.irfftplan * (map_ampspectrum(psgen.psmodel, psgen.noisesignal)[1] .* generate_gaussian_noise(psgen.noisesignal))
end

#
# Compute frequencies in the Fourier Domain
#
@inline function rfftfreq(psgen::PSNoiseGenerator1D)
    return rfftfreq(psgen.noisesignal)
end

#
# Retrieve power-spectrum-scaled fourier noise from signal noise
#
@inline function get_fourier_noise(psgen::PSNoiseGenerator1D, signoise::Vector)::Vector
    return psgen.rfftplan*signoise
end

#
# Retrieve power spectrum from signal noise
#
@inline function get_power_spectrum(psgen::PSNoiseGenerator1D, signoise::Vector)::Vector
    return (abs.(get_fourier_noise(psgen, signoise))).^2

end