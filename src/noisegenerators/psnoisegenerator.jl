export PowerSpectrumNoiseGenerator
export generate_signal_noise
export get_power_spectrum
export get_fourier_noise

"""
    PowerSpectrumNoiseGenerator

This is an abstract data type for power spectrum noise generators.

# Mandatory fields
-`psmodel::AbstractNoisePowerSpectrum`:
-`noisesignal::AbstractNoiseSignal`:
-`rfftplan`:
-`irfftplan`:

# Mandatory methods
-`generate_signal_noise(psgen::PowerSpectrumNoiseGenerator)`: transforms fourier noise to signal noise in position plane for given power spectrum model
-`rfftfreq(psgen::PowerSpectrumNoiseGenerator)`: returns the frequency grid along each dimension in the Fourier plane
-`get_fourier_noise(psgen::PowerSpectrumNoiseGenerator, signoise)`: returns fourier noise for given signal noise
-`get_power_spectrum(psgen::PowerSpectrumNoiseGenerator, signoise:)`: returns power spectrum model retrieved from signal noise
"""

abstract type PowerSpectrumNoiseGenerator end

#
# Convert fourier domain gaussian noise to signal noise with power law and inverse rfft 
#
@inline function generate_signal_noise(psgen::PowerSpectrumNoiseGenerator)::AbstractArray
    return psgen.irfftplan * (map_ampspectrum(psgen.psmodel, psgen.noisesignal)[1] .* generate_gaussian_noise(psgen.noisesignal))
end

#
# Compute frequencies in the Fourier Domain
#
@inline function AbstractFFTs.rfftfreq(psgen::PowerSpectrumNoiseGenerator)
    return rfftfreq(psgen.noisesignal)
end

#
# Retrieve power-spectrum-scaled fourier noise from signal noise
#
@inline function get_fourier_noise(psgen::PowerSpectrumNoiseGenerator, signoise::AbstractArray)::AbstractArray
    return psgen.rfftplan*signoise
end

#
# Retrieve power spectrum from signal noise
#
@inline function get_power_spectrum(psgen::PowerSpectrumNoiseGenerator, signoise::AbstractArray)::AbstractArray
    return (abs.(get_fourier_noise(psgen, signoise))).^2

end