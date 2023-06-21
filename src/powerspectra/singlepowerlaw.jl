export SinglePowerLaw

"""
    SinglePowerLaw <: AbstractPowerSpectrumModel

This is a data type for an N dimensional power law function of the form P(ν)=ν^(index)
"""
struct SinglePowerLaw{N,T} <: AbstractPowerSpectrumModel{N}
    index::T
    function SinglePowerLaw{N}(index::T) where {N,T}
        return new{N,T}(index)
    end
end

@inline fourieranalytic(::SinglePowerLaw) = IsAnalytic()

power_point(model::SinglePowerLaw, ν...) = √(sum(abs2, ν))^model.index

function power_map(model::AbstractPowerSpectrumModel, gridofν::Tuple) 
    if length(gridofν) == 1
        pow = [power_point(model, u) for u in gridofν[1]]
    elseif length(gridofν) == 2
        pow = [power_point(model,u,v) for u in gridofν[1], v in gridofν[2]]
    elseif length(gridofν) == 3
        pow = [power_point(model, u, v, w) for u in gridofν[1], v in gridofν[2], w in gridofν[3]]
    end
    pow[1] = 0
    return pow
end

power_map(model::AbstractPowerSpectrumModel, noisesignal::Union{AbstractNoiseSignal, AbstractContinuousNoiseSignal}) = power_map(model, rfftfreq(noisesignal))
