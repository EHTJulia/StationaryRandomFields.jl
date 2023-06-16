export SinglePowerLaw

"""
    SinglePowerLaw1D <: AbstractNoisePowerSpectrum

This is a data type for an N dimensional power law function of the form P(ν)=amp*ν^(index)
"""

struct SinglePowerLaw <: AbstractNoisePowerSpectrum 
    amp::Float64
    index::Float64
end

