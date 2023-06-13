export AbstractNoisePowerSpectrum
export NoisePowerSpectrum1D
export NoisePowerSpectrum2D


"""
    AbstractNoisePowerSpectrum

This is an abstract data type for noise power spectra.

# Mandatory fields
- `amp::Float64`: amplitude of power spectrum
- `index::Float64`: power law index β 
# Mandatory methods
-`map_ampspectrum_point(psmodel::AbstractNoisePowerSpectrum, singleν::Number)`: maps power spectrum amplitude for a single frequency point
-`map_ampspectrum(psmodel::AbstractNoisePowerSpectrum, gridofν::AbstractArray...)`: maps power spectrum amplitudes for a grid of frequencies
-`map_powerspectrum_point(psmodel::AbstractNoisePowerSpectrum, singleν::Number)`: maps power spectrum for a single frequency point
-`map_powerspectrum(psmodel::AbstractNoisePowerSpectrum, gridofν::AbstractArray...)`: maps power spectrum for a grid of frequencies

"""
abstract type AbstractNoisePowerSpectrum end
Base.ndims(::AbstractNoisePowerSpectrum) = 1

abstract type NoisePowerSpectrum1D <: AbstractNoisePowerSpectrum end
Base.ndims(::NoisePowerSpectrum1D) = 1

abstract type NoisePowerSpectrum2D <: AbstractNoisePowerSpectrum end
Base.ndims(::NoisePowerSpectrum2D) = 2


