"""
    PowerSpectrumNoiseGenerator

This is an abstract data type for power spectrum noise generators.

# Mandatory fields

# Mandatory methods

"""

abstract type PowerSpectrumNoiseGenerator end

@inline function generate_signal_noise(psgen::PowerSpectrumNoiseGenerator)
    return psgen.ifftplan * (map_ampspectrum(psgen.psmodel, psgen.noisesignal) * generate_gaussian_noise(psgen.noisesignal))
end