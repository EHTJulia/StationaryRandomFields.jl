export AbstractNoisePowerSpectrum
export NoisePowerSpectrum1D
export NoisePowerSpectrum2D
export map_ampspectrum
export map_ampspectrum_point
export map_powerspectrum
export map_powerspectrum_point

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


#
# Method to compute given amplitude power law function of single point frequency
#
@inline function map_ampspectrum_point(psmodel::AbstractNoisePowerSpectrum, singleν::Number)::Number #is number data type correct here
    return psmodel.amp * singleν^(psmodel.index/2) 
end

"""
    mapampspectrum
"""

#
# Compute amplitude spectrum of frequency grid
#
@inline function map_ampspectrum(psmodel::AbstractNoisePowerSpectrum, gridofν::AbstractArray...)
    return (gridofν -> psmodel.amp .* gridofν .^ psmodel.index/2).(gridofν)
end

@inline function map_ampspectrum(psmodel::AbstractNoisePowerSpectrum, gridofν::Tuple)
    return (gridofν -> psmodel.amp .* gridofν .^ psmodel.index/2).(gridofν)
end

#
# Compute amplitude spectrum of frequency grid corresponding to signal data
#
@inline function map_ampspectrum(psmodel::AbstractNoisePowerSpectrum, signaldata::AbstractNoiseSignal)
    return map_ampspectrum(psmodel,rfftfreq(signaldata))
end

#
# Compute power law spectrum of single point frequency 
#
@inline function map_powerspectrum_point(psmodel::AbstractNoisePowerSpectrum, singleν::Number)::Number
    return map_ampspectrum_point(psmodel, singleν)^2
end 

#
# Compute amplitude spectrum of frequency grid 
#
@inline function map_powerspectrum(psmodel::AbstractNoisePowerSpectrum, gridofν::AbstractArray...) 
    ampspec = map_ampspectrum(psmodel, gridofν)
    return (ampspec -> ampspec .^ 2).(ampspec)
end

@inline function map_powerspectrum(psmodel::AbstractNoisePowerSpectrum, gridofν::Tuple) 
    ampspec = map_ampspectrum(psmodel, gridofν)
    return (ampspec -> ampspec .^ 2).(ampspec)
end

#
# Compute power spectrum of frequency grid corresponding to signal data
#
@inline function map_powerspectrum(psmodel::AbstractNoisePowerSpectrum, signaldata::NoiseSignal1D)
    return map_powerspectrum(psmodel,rfftfreq(signaldata))
end


