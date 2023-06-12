"""
    AbstractNoisePowerSpectrum

This is an abstract data type for noise power spectra.

# Mandatory fields
- `amp::Float64`: amplitude of power spectrum
- `index::Float64`: power law index β 
# Mandatory methods

"""
abstract type AbstractNoisePowerSpectrum end
Base.ndims(::AbstractNoisePowerSpectrum) = 1

abstract type NoisePowerSpectrum1D <: AbstractNoisePowerSpectrum end
Base.ndims(::NoisePowerSpectrum1D) = 1

abstract type NoisePowerSpectrum2D <: AbstractNoisePowerSpectrum end
Base.ndims(::NoisePowerSpectrum2D) = 2


#
# Method to compute given amplitude power law function of single point signal
#
@inline function map_ampspectrum_point(psmodel::AbstractNoisePowerSpectrum, singleν::Number)::Number #is number data type correct here
    return psmodel.amp * singleν^(psmodel.index) 
end

#
# Method to compute given amplitude power law function of signal grid 
#
@inline function map_ampspectrum(psmodel::AbstractNoisePowerSpectrum, gridofν::Tuple)::Tuple # should i define this for other data types
    return (singleν -> psmodel.amp .* singleν .^ psmodel.index).(singleν)
end

#
# Method to compute given power law function of signal grid 
#
@inline function map_powerspectrum_point(psmodel::AbstractNoisePowerSpectrum, singleν::Number)::Number
    return map_ampspectrum_point(psmodel, singleν)^2

#
# Method to compute given amplitude power law function of signal grid 
#
@inline function map_powerspectrum(psmodel::AbstractNoisePowerSpectrum, gridofν::Tuple)::Tuple # should i define this for other data types
    ampspec = map_ampspectrum(psmodel,gridofν)
    return (ampspec -> ampspec .^ 2).(ampspec)
end
