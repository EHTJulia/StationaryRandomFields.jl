export BrokenPowerLaw
"""
    $(TYPEDEF)

Type for power law with input dimension {N} and profile
```math
    P(\nu) = | \nu | ^{\beta} \end{cases}
```
when inscale < |ν| < outscale. β refers to the input index. 

When |ν| > outscale: P(ν) = 0, and when |ν| < inscale:

```math
    P(\nu) = (inscale) ^{\beta} \end{cases}
```
Can be renormalized, stretched, and rotated via ModelModifier.

# Fields
$(FIELDS)
"""
struct BrokenPowerLaw{N,T} <: AbstractPowerSpectrumModel{N}
    index::T
    inscale::Number
    outscale::Number
    function BrokenPowerLaw{N}(index::T, inscale::Number, outscale::Number) where {N, T}
        if outscale <= inscale
            throw(ArgumentError("Value of 'outscale' must be greater than 'inscale.'"))
        end
        new{N, typeof(index)}(index, inscale, outscale)
    end
end

function power_point(model::BrokenPowerLaw, ν...)
    if √sum(abs2, ν) <= model.inscale
        return model.inscale^model.index
    elseif √sum(abs2, ν) >= model.outscale
        return 0.
    else
        return √(sum(abs2, ν))^model.index
    end
end



    
