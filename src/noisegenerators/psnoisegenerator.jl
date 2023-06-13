export PowerSpectrumNoiseGenerator
export generate_signal_noise

"""
    PowerSpectrumNoiseGenerator

This is an abstract data type for power spectrum noise generators.

# Mandatory fields

# Mandatory methods
- `generate_signal_noise(psgen::PowerSpectrumNoiseGenerator)`: transforms fourier noise to signal noise for given power spectrum model

"""

abstract type PowerSpectrumNoiseGenerator end

