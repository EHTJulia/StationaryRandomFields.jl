"""A periodic Cartesian grid; spacing is measured in signal-coordinate units."""
struct FieldGrid{N,T<:AbstractFloat}
    dims::NTuple{N,Int}
    spacing::NTuple{N,T}
end
function FieldGrid(dims::NTuple{N,Integer}; spacing=ntuple(_ -> 1.0, N)) where {N}
    N > 0 && all(>(0), dims) || throw(ArgumentError("Invalid grid dimensions."))
    length(spacing) == N || throw(DimensionMismatch("Grid spacing has the wrong rank."))
    s = float.(promote(spacing...))
    all(v -> isfinite(v) && v > 0, s) || throw(ArgumentError("Invalid grid spacing."))
    return FieldGrid(Int.(dims), s)
end
Base.size(g::FieldGrid) = g.dims
Base.ndims(::FieldGrid{N}) where {N} = N

"""A nonnegative power spectrum with `N` frequency coordinates."""
abstract type AbstractPowerSpectrumModel{N} end
Base.ndims(::AbstractPowerSpectrumModel{N}) where {N} = N

"""A constant discrete power spectrum."""
struct ConstantPowerSpectrum{N,T} <: AbstractPowerSpectrumModel{N}
    power::T
    function ConstantPowerSpectrum{N}(p::T; check_args::Bool=true) where {N,T}
        (!check_args || (isfinite(p) && p >= 0)) ||
            throw(ArgumentError("Power must be finite and nonnegative."))
        N > 0 || throw(ArgumentError("Spectrum rank must be positive."))
        return new{N,T}(p)
    end
end

"""A radial power law. The DC power is specified separately and defaults to zero."""
struct SinglePowerLaw{N,T,S} <: AbstractPowerSpectrumModel{N}
    index::T
    dc::S
    function SinglePowerLaw{N}(
        index::T; dc::S=zero(index), check_args::Bool=true
    ) where {N,T,S}
        N > 0 && (!check_args || isfinite(index)) ||
            throw(ArgumentError("Invalid power law."))
        (!check_args || (isfinite(dc) && dc >= 0)) ||
            throw(ArgumentError("Invalid DC power."))
        return new{N,T,S}(index, dc)
    end
end

"""A power law flat below `inscale` and zero at or above `outscale`."""
struct SaturatedPowerLaw{N,T} <: AbstractPowerSpectrumModel{N}
    index::T
    inscale::T
    outscale::T
end
function SaturatedPowerLaw{N}(index, inscale, outscale; check_args::Bool=true) where {N}
    β, lo, hi = float.(promote(index, inscale, outscale))
    N > 0 && (!check_args || (all(isfinite, (β, lo, hi)) && 0 < lo < hi)) ||
        throw(ArgumentError("Require finite index and 0 < inscale < outscale."))
    return SaturatedPowerLaw{N,typeof(β)}(β, lo, hi)
end

function check_frequency(p, ν)
    length(ν) == ndims(p) || throw(DimensionMismatch("Wrong frequency rank."))
    return nothing
end
function power_point(p::ConstantPowerSpectrum, ν...)
    check_frequency(p, ν)
    return p.power
end
function power_point(p::SinglePowerLaw, ν...)
    check_frequency(p, ν)
    r2 = sum(abs2, ν)
    # Keep the unused DC branch finite during automatic differentiation.
    r = sqrt(ifelse(iszero(r2), one(r2), r2))
    return ifelse(iszero(r2), p.dc, r^p.index)
end
function power_point(p::SaturatedPowerLaw, ν...)
    check_frequency(p, ν)
    r2 = sum(abs2, ν)
    r = sqrt(max(r2, p.inscale^2))
    return ifelse(r2 >= p.outscale^2, zero(r), r^p.index)
end
# Zero-power modes stay fixed at zero; use a zero derivative at this boundary.
function power_amplitude(power)
    safe = ifelse(iszero(power), one(power), power)
    return ifelse(iszero(power), zero(power), sqrt(safe))
end
amplitude_point(p::AbstractPowerSpectrumModel, ν...) = power_amplitude(power_point(p, ν...))

struct SumPowerSpectrum{N,P<:Tuple} <: AbstractPowerSpectrumModel{N}
    models::P
end
function Base.:+(
    a::AbstractPowerSpectrumModel{N}, b::AbstractPowerSpectrumModel{N}
) where {N}
    return SumPowerSpectrum{N,typeof((a, b))}((a, b))
end
power_point(p::SumPowerSpectrum, ν...) = sum(m -> power_point(m, ν...), p.models)

abstract type ModelModifier end
struct Renormalize{T<:Real} <: ModelModifier
    scale::T
    function Renormalize(s::T) where {T<:Real}
        isfinite(s) && s >= 0 || throw(ArgumentError("Power scaling must be nonnegative."))
        return new{T}(s)
    end
end
struct Stretch{N,T<:Real} <: ModelModifier
    factors::NTuple{N,T}
end
function Stretch(s::Real...)
    f = promote(s...)
    !isempty(f) && all(v -> isfinite(v) && v > 0, f) ||
        throw(ArgumentError("Stretch factors must be positive."))
    return Stretch(f)
end
struct Rotate{N,T} <: ModelModifier
    rows::NTuple{N,NTuple{N,T}}
end
function Rotate(angle::Real)
    s, c = sincos(angle)
    return Rotate(((c, -s), (s, c)))
end
function Rotate(q::AbstractMatrix{<:Real})
    n, m = size(q)
    n == m && n > 0 || throw(DimensionMismatch("Rotation must be square."))
    isapprox(q' * q, I; atol=100eps(float(eltype(q))), rtol=100eps(float(eltype(q)))) &&
    isapprox(det(q), one(eltype(q)); atol=100eps(float(eltype(q)))) ||
        throw(ArgumentError("Rotation must be a proper orthogonal matrix."))
    return Rotate(ntuple(i -> ntuple(j -> q[i, j], n), n))
end
struct ModifiedPowerSpectrumModel{N,P<:AbstractPowerSpectrumModel{N},M<:Tuple} <:
       AbstractPowerSpectrumModel{N}
    model::P
    modifiers::M
end
function modify(p::AbstractPowerSpectrumModel{N}, ms::ModelModifier...) where {N}
    for m in ms
        m isa Stretch &&
            length(m.factors) != N &&
            throw(DimensionMismatch("Stretch has the wrong rank."))
        m isa Rotate &&
            length(m.rows) != N &&
            throw(DimensionMismatch("Rotation has the wrong rank."))
    end
    return ModifiedPowerSpectrumModel{N,typeof(p),typeof(ms)}(p, ms)
end
function power_point(p::ModifiedPowerSpectrumModel, ν...)
    check_frequency(p, ν)
    return modified_power(p.model, p.modifiers, ν)
end
modified_power(p, ::Tuple{}, ν) = power_point(p, ν...)
function modified_power(p, ms::Tuple, ν)
    return modified_power(p, Base.front(ms), last(ms), ν)
end
modified_power(p, ms, m::Renormalize, ν) = m.scale * modified_power(p, ms, ν)
modified_power(p, ms, m::Stretch, ν) = modified_power(p, ms, ν ./ m.factors)
function modified_power(p, ms, m::Rotate, ν)
    modified_power(p, ms, map(row -> sum(row .* ν), m.rows))
end
renormed(p::AbstractPowerSpectrumModel, s::Real) = modify(p, Renormalize(s))
stretched(p::AbstractPowerSpectrumModel, s::Real...) = modify(p, Stretch(s...))
rotated(p::AbstractPowerSpectrumModel, q) = modify(p, Rotate(q))
Base.:*(p::AbstractPowerSpectrumModel, s::Real) = renormed(p, s)
Base.:*(s::Real, p::AbstractPowerSpectrumModel) = renormed(p, s)
function Base.:/(p::AbstractPowerSpectrumModel, s::Real)
    isfinite(s) && s > 0 || throw(ArgumentError("Power divisor must be positive."))
    return renormed(p, inv(s))
end

function frequency_axes(g::FieldGrid{N}, (::Type{T})=Float64) where {N,T}
    return ntuple(N) do d
        n = g.dims[d]
        [T(k <= (n - 1) ÷ 2 ? k : k - n) / T(n * g.spacing[d]) for k in 0:(n - 1)]
    end
end

# Copy metadata to the destination backend once, before a compiled transform.
function frequency_coordinates(g::FieldGrid{N}, prototype; reflected=false) where {N}
    axes = frequency_axes(g, eltype(prototype))
    return ntuple(N) do d
        host = reshape(
            (reflected ? axes[d][[1; collect(length(axes[d]):-1:2)]] : axes[d]),
            ntuple(i -> i == d ? g.dims[d] : 1, N),
        )
        device = similar(prototype, eltype(prototype), size(host))
        copyto!(device, host)
    end
end
function power_map!(out, p::AbstractPowerSpectrumModel, coordinates::Tuple)
    length(coordinates) == ndims(p) || throw(DimensionMismatch("Wrong grid rank."))
    out .= power_point.(Ref(p), coordinates...)
    return out
end
function power_map(p::AbstractPowerSpectrumModel, axes::Tuple)
    length(axes) == ndims(p) || throw(DimensionMismatch("Wrong grid rank."))
    coords = ntuple(
        d -> reshape(axes[d], ntuple(i -> i == d ? length(axes[d]) : 1, length(axes))),
        length(axes),
    )
    return power_point.(Ref(p), coords...)
end
power_map(p::AbstractPowerSpectrumModel, g::FieldGrid) = power_map(p, frequency_axes(g))
amplitude_map(p::AbstractPowerSpectrumModel, grid) = sqrt.(power_map(p, grid))
function amplitude_map!(out, p::AbstractPowerSpectrumModel, coordinates::Tuple)
    power_map!(out, p, coordinates)
    out .= sqrt.(out)
    return out
end
