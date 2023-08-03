export PowerSpectrumNoiseGenerator
export generate_signal_noise
export get_power_spectrum
export get_fourier_noise

"""
    PowerSpectrumNoiseGenerator

This is an abstract data type for power spectrum noise generators.

**Mandatory methods**
- [`generate_signal_noise`](@ref): Transforms fourier noise to signal noise in position plane for given power spectrum model
- [`get_fourier_noise`](@ref): Computes fourier noise for given signal domain noise
- [`get_power_spectrum`](@ref): Returns power spectrum model retrieved from given signal domain noise
"""

abstract type PowerSpectrumNoiseGenerator end

"""
    generate_signal_noise(psgen::PowerSpectrumNoiseGenerator)

Converts fourier domain gaussian noise to signal noise, corresponding to the continuous noise signal and power spectrum
contained in psgen.
Scaling the gaussian fourier noise by the power law function and then transforming to the signal domain with an inverse 
rfft generates the signal noise.
"""

@inline function generate_signal_noise(psgen::PowerSpectrumNoiseGenerator)::AbstractArray
    return psgen.noisesignal.invplan * (amplitude_map(psgen.psmodel, psgen.noisesignal) .* generate_gaussian_noise(psgen.noisesignal))
end

"""
    get_fourier_noise(psgen::PowerSpectrumNoiseGenerator, signoise::AbstractArray)

Retrieve power-spectrum-scaled fourier noise from given signal noise by transforming to fourier space.
"""
@inline function get_fourier_noise(psgen::PowerSpectrumNoiseGenerator, signoise::AbstractArray)::AbstractArray
    return psgen.noisesignal.plan*signoise
end

"""
    get_power_spectrum(psgen::PowerSpectrumNoiseGenerator, signoise::AbstractArray)

Retrieve power spectrum from given signal noise.
"""
@inline function get_power_spectrum(psgen::PowerSpectrumNoiseGenerator, signoise::AbstractArray)::AbstractArray
    return .5 * (abs.(get_fourier_noise(psgen, signoise))).^2
end

@inline function AbstractFFTs.rfftfreq(psgen::PowerSpectrumNoiseGenerator)
    return rfftfreq(psgen.noisesignal)
end