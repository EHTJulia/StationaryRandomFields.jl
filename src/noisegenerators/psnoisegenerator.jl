export PowerSpectrumNoiseGenerator
export generate_signal_noise

"""
    PowerSpectrumNoiseGenerator

This is an abstract data type for power spectrum noise generators.

# Mandatory fields

# Mandatory methods

"""

abstract type PowerSpectrumNoiseGenerator end

@inline function generate_signal_noise(psgen::PowerSpectrumNoiseGenerator)
    ampspec = map_ampspectrum(psgen.psmodel, psgen.noisesignal)
    ampsec[1] = 0
    return psgen.ifftplan * (ampspec * generate_gaussian_noise(psgen.noisesignal))
end