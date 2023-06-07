"""
    AbstractNoisePowerSpectrum

This is an abstract data type for noise power spectra.

# Mandatory fields

# Mandatory methods

"""
abstract type AbstractNoisePowerSpectrum end
Base.ndims(::AbstractNoisePowerSpectrum) = 1

abstract type NoisePowerSpectrum1D <: AbstractNoisePowerSpectrum end
Base.ndims(::NoisePowerSpectrum1D) = 1

abstract type NoisePowerSpectrum2D <: AbstractNoisePowerSpectrum end
Base.ndims(::NoisePowerSpectrum1D) = 2