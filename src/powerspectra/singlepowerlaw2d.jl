export SinglePowerLaw2D
export map_ampspectrum
export map_ampspectrum_point
export map_powerspectrum
export map_powerspectrum_point

"""
    SinglePowerLaw2D <: NoisePowerSpectrum2D

This is a data type for a 2 dimensional power law function of the form P(ν)=amp*|ν|^(index)
"""

struct SinglePowerLaw2D <: NoisePowerSpectrum2D 
    amp::Float64
    index::Float64
end


#
# Method to compute given amplitude power law function of single 2D point frequency
#
@inline function map_ampspectrum_point(psmodel::NoisePowerSpectrum2D, singleν::Tuple)::Number 
    return psmodel.amp * (singleν[1]^2. + singleν[2]^2.)^(-psmodel.index/2)
end

"""
    mapampspectrum
"""

#
# Compute amplitude spectrum of 2D frequency grid
#
@inline function map_ampspectrum(psmodel::NoisePowerSpectrum2D, normgrid::AbstractArray...)
    ampspec = (normgrid -> psmodel.amp .* normgrid .^ (-psmodel.index/2)).(normgrid)
    ampspec[1][1] = 0
    return ampspec
end

@inline function map_ampspectrum(psmodel::NoisePowerSpectrum2D, gridofν::Tuple)
    return map_ampspectrum(psmodel, freq_norm(gridofν))
end

#
# Compute amplitude spectrum of frequency grid corresponding to 2D signal data
#
@inline function map_ampspectrum(psmodel::NoisePowerSpectrum2D, signaldata::NoiseSignal2D)
    return map_ampspectrum(psmodel, freq_norm(signaldata))
end


#
# Compute power law spectrum of single 2D point frequency 
#
@inline function map_powerspectrum_point(psmodel::NoisePowerSpectrum2D, singleν::Tuple)::Number
    return map_ampspectrum_point(psmodel, singleν)^2
end 

"""
    map_powerspectrum
"""

#
# Compute amplitude spectrum of 2D frequency grid 
#
@inline function map_powerspectrum(psmodel::NoisePowerSpectrum2D, data)
    ampspec = map_ampspectrum(psmodel, data)
    return (ampspec -> ampspec .^ 2).(ampspec)
end
