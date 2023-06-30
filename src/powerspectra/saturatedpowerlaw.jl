export SaturatedPowerLaw
@doc raw"""
    $(TYPEDEF)

Type for power law with input dimension {N} and profile: ``P(\nu) = |\nu| ^{\beta} \end{cases}``
when inscale < |ν| < outscale.
When |ν| > outscale: P(ν) = 0.
When |ν| < inscale: ``P(\\nu) = (inscale) ^{\\beta} \\end{cases}.``

β refers to the input index. Can be renormalized, stretched, and rotated via ModelModifier.

# Fields
$(FIELDS)
"""
struct SaturatedPowerLaw{N,T} <: AbstractPowerSpectrumModel{N}
    index::T
    inscale::Number
    outscale::Number
    function SaturatedPowerLaw{N}(index::T, inscale::Number, outscale::Number) where {N,T}
        if outscale <= inscale
            throw(ArgumentError("Value of 'outscale' must be greater than 'inscale.'"))
        end
        new{N,typeof(index)}(index, inscale, outscale)
    end
end

@inline fourieranalytic(::SinglePowerLaw) = IsAnalytic()

function power_point(model::SaturatedPowerLaw, ν...)
    @assert length(ν) == ndims(model)

    norm = √sum(abs2, ν)
    if norm <= model.inscale
        return model.inscale^model.index
    elseif norm >= model.outscale
        return 0.0
    else
        return norm^model.index
    end
end
