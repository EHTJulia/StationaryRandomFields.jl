export SinglePowerLaw2D

"""
    SinglePowerLaw2D <: NoisePowerSpectrum1D

This is a data type for a 2 dimensional power law function of the form P(ν)=amp*|ν|^(index)
"""

struct SinglePowerLaw2D <: NoisePowerSpectrum2D 
    amp::Float64
    index::Float64
end
