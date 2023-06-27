export SinglePowerLaw
"""
    $(TYPEDEF)

Type for single power law with input dimension {N} and profile
```math
    P(\nu) = | \nu | ^{\beta} \end{cases}
```
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

power_point(model::SinglePowerLaw, ν...) = √(sum(abs2, ν))^model.index

