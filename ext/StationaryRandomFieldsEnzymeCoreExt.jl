module StationaryRandomFieldsEnzymeCoreExt

using StationaryRandomFields: apply, PreparedWavelet, wavelet_analysis
using EnzymeCore: Const, Duplicated
using EnzymeCore.EnzymeRules:
    EnzymeRules, FwdConfigWidth, RevConfigWidth, needs_primal, needs_shadow, AugmentedReturn

function check_active_buffer(config, x)
    if EnzymeRules.runtime_activity(config) && x.val === x.dval
        throw(ArgumentError("Wavelet rules require an independent derivative buffer."))
    end
    return nothing
end

# Wavelet filters and packing are fixed metadata. Orthogonal analysis is the
# transpose of synthesis; no signal-space density or inverse Jacobian is used.
function EnzymeRules.forward(
    config::FwdConfigWidth{1},
    ::Const{typeof(apply)},
    RT::Type,
    p::Const{<:PreparedWavelet},
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
    p::Const{<:PreparedWavelet},
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
    p::Const{<:PreparedWavelet},
    x::Duplicated{<:Array},
)
    x.dval .+= wavelet_analysis(p.val, dy)
    fill!(dy, zero(eltype(dy)))
    return (nothing, nothing)
end

end
