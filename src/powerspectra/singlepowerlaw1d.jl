export SinglePowerLaw1D

"""
    SinglePowerLaw1D <: NoisePowerSpectrum1D

This is a data type for a 1 dimensional power law function of the form P(ν)=amp*ν^(index)
"""

struct SinglePowerLaw1D <: NoisePowerSpectrum1D 
    amp::Float64
    index::Float64
end




