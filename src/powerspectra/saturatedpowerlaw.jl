export SaturatedPowerLaw
"""
    $(TYPEDEF)

Type for power law with input dimension {N} and profile
```math
    P(\nu) = | \nu | ^{\beta} \end{cases}.
```
when inscale < |ν| < outscale. β refers to the input index. 

When |ν| > outscale: P(ν) = 0. 
When |ν| < inscale:

```math
    P(\nu) = (inscale) ^{\beta} \end{cases}.
```
Can be renormalized, stretched, and rotated via ModelModifier.

# Fields
$(FIELDS)
"""
struct SaturatedPowerLaw{N,T} <: AbstractPowerSpectrumModel{N}
    index::T
    inscale::Number
    outscale::Number
    function SaturatedPowerLaw{N}(index::T, inscale::Number, outscale::Number) where {N, T}
        if outscale <= inscale
            throw(ArgumentError("Value of 'outscale' must be greater than 'inscale.'"))
        end
        new{N, typeof(index)}(index, inscale, outscale)
    end
end

@inline fourieranalytic(::SinglePowerLaw) = IsAnalytic()

function power_point(model::SaturatedPowerLaw, ν...)
    norm = √sum(abs2, ν)
    if norm <= model.inscale
        return model.inscale^model.index
    elseif norm >= model.outscale
        return 0.
    else
        return norm^model.index
    end
end



    
