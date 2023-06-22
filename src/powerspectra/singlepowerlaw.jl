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
    prod = collect(Iterators.product(gridofν...))
    pow = zeros(size(prod)...)
    for i in eachindex(prod)[2:end]
        pow[i] = power_point(model,prod[i]...)
    end
    return pow
end

power_map(model::AbstractPowerSpectrumModel, noisesignal::Union{AbstractNoiseSignal, AbstractContinuousNoiseSignal}) = power_map(model, rfftfreq(noisesignal))
