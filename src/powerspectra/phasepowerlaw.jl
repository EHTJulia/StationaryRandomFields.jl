export PhaseScreenPowerLaw
@doc raw"""
    $(TYPEDEF)

D is observer screen distance

# Fields
$(FIELDS)
"""
struct PhaseScreenPowerLaw{N,T, F<:Function} <: AbstractPowerSpectrumModel{2}
    α::T
    P_ϕ::F
    r_in::Number
    Qbar::Number
    D::Number
    Vx_km_per_s::Number
    Vy_km_per_s::Number

    function PhaseScreenPowerLaw(index::T, P_ϕ::F, r_in::Number, Qbar::Number, D = 2.82 * 3.086e21, Vx_km_per_s=50.0, Vy_km_per_s=0.0) where {T, F}
        new{2,typeof(index), typeof(P_ϕ)}(index, P_ϕ, r_in, Qbar, D, Vx_km_per_s, Vy_km_per_s)
    end
end

@inline fourieranalytic(::PhaseScreenPowerLaw) = IsAnalytic()

function power_point(model::PhaseScreenPowerLaw, ν...)
    @assert length(ν) == 2

    q = √sum(abs2, ν) + 1e-12/model.r_in 
    ϕ = atan(ν[2], ν[1])
    return model.Qbar * (q*model.r_in)^(-(model.α + 2.0)) * exp(-(q * model.r_in)^2) * model.P_ϕ(ϕ)
end

function amplitude_point(model::PhaseScreenPowerLaw, t_hr, N, FOV, dq, ν...)
    x_offset_pixels = (model.Vx_km_per_s * 1.e5) * (t_hr*3600.) / (FOV/float(N))
    y_offset_pixels = (model.Vy_km_per_s * 1.e5) * (t_hr*3600.) / (FOV/float(N))
    s = ν[1]
    t = ν[2]

    return √(power_point(model, dq*s, dq*t)) * exp(2.0*π*1im*(s*x_offset_pixels + t*y_offset_pixels)/float(N))
end

function amplitude_map(model::PhaseScreenPowerLaw, noisesignal::Union{AbstractNoiseSignal,AbstractContinuousNoiseSignal}, t_hr = 0.) 
    N = noisesignal.dims[1]
    FOV = N * model.D
    dq = 2*π/FOV
    gridofν = rfftfreq(noisesignal)

    prod = collect(Iterators.product(gridofν...))
    amp = zeros(size(prod)...)
    for i in eachindex(prod)[2:end]
        amp[i] = amplitude_point(model, t_hr, N, FOV, dq, prod[i]...)
    end
    return amp
end
