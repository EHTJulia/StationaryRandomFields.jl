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

#
# Method to compute given amplitude power law function of single 2D point frequency
#
@inline function map_ampspectrum_point(psmodel::AbstractNoisePowerSpectrum, singleν::Tuple)::Number 
    return psmodel.amp * (singleν[1]^2. + singleν[2]^2.)^(-psmodel.index/2)
end

@inline function map_ampspectrum_point(psmodel::AbstractNoisePowerSpectrum, singleν::Number)::Number 
    return psmodel.amp * singleν^(-psmodel.index/2) 
end

"""
    mapampspectrum
"""

#
# Compute amplitude spectrum of 2D frequency grid
#
@inline function map_ampspectrum(psmodel::AbstractNoisePowerSpectrum, normgrid::AbstractArray...)
    ampspec = (normgrid -> psmodel.amp .* normgrid .^ (-psmodel.index/2)).(normgrid)
    ampspec[1][1] = 0
    return ampspec
end

@inline function map_ampspectrum(psmodel::AbstractNoisePowerSpectrum, gridofν::Tuple)
    return map_ampspectrum(psmodel, freq_norm(gridofν))
end

#
# Compute amplitude spectrum of frequency grid corresponding to 2D signal data
#
@inline function map_ampspectrum(psmodel::AbstractNoisePowerSpectrum, signaldata::AbstractNoiseSignal)
    return map_ampspectrum(psmodel, freq_norm(signaldata))
end



#
# Compute power law spectrum of single 2D point frequency 
#
@inline function map_powerspectrum_point(psmodel::AbstractNoisePowerSpectrum, singleν)::Number
    return map_ampspectrum_point(psmodel, singleν)^2
end 

#
# Compute amplitude spectrum of 2D frequency grid 
#
@inline function map_powerspectrum(psmodel::AbstractNoisePowerSpectrum, data)
    ampspec = map_ampspectrum(psmodel, data)
    return (ampspec -> ampspec .^ 2).(ampspec)
end
