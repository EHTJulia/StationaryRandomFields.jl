module StationaryRandomFieldsEnzymeCoreExt

using StationaryRandomFields: apply, PreparedWavelet, PreparedHartley, wavelet_analysis

const PreparedLinearSynthesis = Union{PreparedWavelet, PreparedHartley}
transpose_apply(p::PreparedWavelet, x) = wavelet_analysis(p, x)
transpose_apply(p::PreparedHartley, x) = apply(p, x)
using EnzymeCore: Const, Duplicated
using EnzymeCore.EnzymeRules:
    EnzymeRules, FwdConfigWidth, RevConfigWidth, needs_primal, needs_shadow, AugmentedReturn

function check_active_buffer(config, x)
    if EnzymeRules.runtime_activity(config) && x.val === x.dval
        throw(
            ArgumentError(
                "Linear synthesis rules require an independent derivative buffer."
            ),
        )
    end
    return nothing
end

# Plans, wavelet filters, and packing are fixed metadata. The reverse rule
# applies the linear transpose without introducing a density or Jacobian term.
function EnzymeRules.forward(
        config::FwdConfigWidth{1},
        ::Const{typeof(apply)},
        RT::Type,
        p::Const{<:PreparedLinearSynthesis},
        x::Duplicated{<:Array},
    )
    check_active_buffer(config, x)
    primal = needs_primal(config) ? apply(p.val, x.val) : nothing
    tangent = needs_shadow(config) ? apply(p.val, x.dval) : nothing
    if needs_primal(config) && needs_shadow(config)
        return Duplicated(primal, tangent)
    end
    return needs_primal(config) ? primal : tangent
end
function EnzymeRules.augmented_primal(
        config::RevConfigWidth{1},
        ::Const{typeof(apply)},
        RT::Type,
        p::Const{<:PreparedLinearSynthesis},
        x::Duplicated{<:Array},
    )
    check_active_buffer(config, x)
    y = needs_primal(config) ? apply(p.val, x.val) : nothing
    dy = zero(x.val)
    return AugmentedReturn(
        needs_primal(config) ? y : nothing, needs_shadow(config) ? dy : nothing, dy
    )
end
function EnzymeRules.reverse(
        config::RevConfigWidth{1},
        ::Const{typeof(apply)},
        RT::Type,
        dy,
        p::Const{<:PreparedLinearSynthesis},
        x::Duplicated{<:Array},
    )
    x.dval .+= transpose_apply(p.val, dy)
    fill!(dy, zero(eltype(dy)))
    return (nothing, nothing)
end

# A constant first transform can still need a shadow when its output is combined
# with active spectrum parameters later in a covariance action.
function EnzymeRules.forward(
        config::FwdConfigWidth{1},
        ::Const{typeof(apply)},
        RT::Type,
        p::Const{<:PreparedLinearSynthesis},
        x::Const{<:Array},
    )
    primal = needs_primal(config) ? apply(p.val, x.val) : nothing
    tangent = needs_shadow(config) ? zero(x.val) : nothing
    if needs_primal(config) && needs_shadow(config)
        return Duplicated(primal, tangent)
    end
    return needs_primal(config) ? primal : tangent
end
function EnzymeRules.augmented_primal(
        config::RevConfigWidth{1},
        ::Const{typeof(apply)},
        RT::Type,
        p::Const{<:PreparedLinearSynthesis},
        x::Const{<:Array},
    )
    primal = needs_primal(config) ? apply(p.val, x.val) : nothing
    shadow = needs_shadow(config) ? zero(x.val) : nothing
    return AugmentedReturn(primal, shadow, nothing)
end
function EnzymeRules.reverse(
        config::RevConfigWidth{1},
        ::Const{typeof(apply)},
        RT::Type,
        tape,
        p::Const{<:PreparedLinearSynthesis},
        x::Const{<:Array},
    )
    return (nothing, nothing)
end

end
