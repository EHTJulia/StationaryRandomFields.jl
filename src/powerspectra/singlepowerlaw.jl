export SinglePowerLaw

"""
    SinglePowerLaw <: AbstractNoisePowerSpectrum

This is a data type for an N dimensional power law function of the form P(ν)=ν^(index)
"""
struct SinglePowerLaw{N,T} <: AbstractPowerSpectrumModel{N}
    index::T
    function SinglePowerLaw{N}(index::T) where {N,T}
        return new{N,T}(index)
    end
end

@inline fourieranalytic(::SinglePowerLaw) = IsAnalytic()

power_point(model::SinglePowerLaw, ν...) = sqrt(sum(abs2, ν))^model.index
