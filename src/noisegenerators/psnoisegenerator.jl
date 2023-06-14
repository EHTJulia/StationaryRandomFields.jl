export PowerSpectrumNoiseGenerator

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

