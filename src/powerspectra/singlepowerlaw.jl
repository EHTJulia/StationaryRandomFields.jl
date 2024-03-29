export SinglePowerLaw
@doc raw"""
    $(TYPEDEF)

Type for single power law with input dimension {N} and profile ``P(\nu) = |\nu| ^{\beta} \end{cases}``
where β is the input index and ν is frequency. Can be renormalized, stretched, and rotated via ModelModifier.

# Fields
$(FIELDS)
"""
struct SinglePowerLaw{N,T} <: AbstractPowerSpectrumModel{N}
    index::T
    function SinglePowerLaw{N}(index::T) where {N,T}
        return new{N,T}(index)
    end
end

@inline fourieranalytic(::SinglePowerLaw) = IsAnalytic()

function power_point(model::SinglePowerLaw, ν...)
    @assert length(ν) == ndims(model)

    return √(sum(abs2, ν))^model.index
end
