export SinglePowerLaw1D
export map_ampspectrum
export map_ampspectrum_point
export map_powerspectrum
export map_powerspectrum_point

"""
    SinglePowerLaw1D <: NoisePowerSpectrum1D

This is a data type for a 1 dimensional power law function of the form P(ν)=amp*ν^(index)
"""

struct SinglePowerLaw1D <: NoisePowerSpectrum1D 
    amp::Float64
    index::Float64
end


#
# Method to compute given amplitude power law function of single point frequency
#
@inline function map_ampspectrum_point(psmodel::NoisePowerSpectrum1D, singleν::Number)::Number 
    return psmodel.amp * singleν^(-psmodel.index/2) 
end

"""
    mapampspectrum
"""

#
# Compute amplitude spectrum of frequency grid
#
@inline function map_ampspectrum(psmodel::NoisePowerSpectrum1D, gridofν::AbstractArray...)
    ampspec = (gridofν -> psmodel.amp .* gridofν .^ (-psmodel.index/2)).(gridofν)
    ampspec[1][1] = 0
    return ampspec
end

@inline function map_ampspectrum(psmodel::NoisePowerSpectrum1D, gridofν::Tuple)
    ampspec = (gridofν -> psmodel.amp .* gridofν .^ (-psmodel.index/2)).(gridofν)
    ampspec[1][1] = 0
    return ampspec
end

#
# Compute amplitude spectrum of frequency grid corresponding to signal data
#
@inline function map_ampspectrum(psmodel::NoisePowerSpectrum1D, signaldata::NoiseSignal1D)
    ampspec = map_ampspectrum(psmodel,rfftfreq(signaldata))
    ampspec[1][1]=0
    return ampspec
end

#
# Compute power law spectrum of single point frequency 
#
@inline function map_powerspectrum_point(psmodel::NoisePowerSpectrum1D, singleν::Number)::Number
    return map_ampspectrum_point(psmodel, singleν)^2
end 

"""
    map_powerspectrum
"""

#
# Compute amplitude spectrum of frequency grid 
#
@inline function map_powerspectrum(psmodel::NoisePowerSpectrum1D, gridofν::AbstractArray...) 
    ampspec = map_ampspectrum(psmodel, gridofν)
    return (ampspec -> ampspec .^ 2).(ampspec)
end

@inline function map_powerspectrum(psmodel::NoisePowerSpectrum1D, gridofν::Tuple) 
    ampspec = map_ampspectrum(psmodel, gridofν)
    return (ampspec -> ampspec .^ 2).(ampspec)
end

#
# Compute power spectrum of frequency grid corresponding to signal data
#
@inline function map_powerspectrum(psmodel::NoisePowerSpectrum1D, signaldata::NoiseSignal1D)
    return map_powerspectrum(psmodel,rfftfreq(signaldata))
end



