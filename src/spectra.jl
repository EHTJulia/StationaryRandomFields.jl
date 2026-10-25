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
        # Homogeneous parameter storage preserves mixed-precision Enzyme derivatives.
        β, dcvalue = promote(index, dc)
        return new{N,typeof(β),typeof(dcvalue)}(β, dcvalue)
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

"""
    MaternPS{N}(rho, nu; check_args = true)

A Matérn spectral shape with positive range `rho` and smoothness `nu`.
Frequencies are in cycles per coordinate unit and `kappa = sqrt(8nu)/rho`.
Power is `(1 + (2pi*f/kappa)^2)^(-nu-N/2)`, with radial frequency `f`.
DC power is one; this is not a unit-variance or continuum-normalized density.
Use `renormed` to set a power scale. The finite grid represents a periodic field.
"""
struct MaternPS{N,T} <: AbstractPowerSpectrumModel{N}
    rho::T
    nu::T
    function MaternPS{N}(rho::Real, nu::Real; check_args::Bool=true) where {N}
        r, v = float.(promote(rho, nu))
        N > 0 || throw(ArgumentError("Spectrum rank must be positive."))
        (!check_args || (isfinite(r) && r > 0 && isfinite(v) && v > 0)) ||
            throw(ArgumentError("Matérn range and smoothness must be finite and positive."))
        return new{N,typeof(r)}(r, v)
    end
end

"""
    MarkovPS{N}(scales::Tuple; check_args = true)

A Markov-type spectral shape `1 / (1 + sum((rho[j]*k)^(2j)))`, where
`k = 2pi*sqrt(sum(abs2, frequencies))`. Scales are in coordinate units,
finite and nonnegative, with at least one positive entry. DC power is one.
A single scale gives a Lorentzian spectrum. This samples a continuum spectral
shape on a periodic grid; it does not imply a sparse lattice precision matrix.
The model does not normalize marginal variance. Use `renormed` for power scaling.
"""
struct MarkovPS{N,M,T} <: AbstractPowerSpectrumModel{N}
    scales::NTuple{M,T}
    function MarkovPS{N}(scales::Tuple; check_args::Bool=true) where {N}
        N > 0 && !isempty(scales) || throw(ArgumentError("Rank and order must be positive."))
        all(x -> x isa Real, scales) || throw(ArgumentError("Scales must be real."))
        s = float.(promote(scales...))
        (!check_args || (all(x -> isfinite(x) && x >= 0, s) && any(>(0), s))) ||
            throw(ArgumentError("Require finite nonnegative scales with a positive entry."))
        return new{N,length(s),eltype(s)}(s)
    end
end

function power_point(p::MaternPS{N}, f::Vararg{Any,K}) where {N,K}
    check_frequency(p, f)
    r2 = sum(abs2, f)
    angular = oftype(p.rho, 2pi)
    return (one(r2) + (angular * p.rho)^2 * r2 / (8 * p.nu))^(-p.nu - oftype(p.nu, N) / 2)
end
function power_point(p::MarkovPS, f::Vararg{Any,K}) where {K}
    check_frequency(p, f)
    r2 = sum(abs2, f)
    k2 = oftype(first(p.scales), 2pi)^2 * r2
    terms = ntuple(j -> (p.scales[j]^2 * k2)^j, length(p.scales))
    return inv(one(k2) + sum(terms))
end

parameter_type(::MaternPS{N,T}) where {N,T} = T
parameter_type(::MarkovPS{N,M,T}) where {N,M,T} = T
function convert_parameters(::Type{T}, p::MaternPS{N}) where {T,N}
    return MaternPS{N}(T(p.rho), T(p.nu); check_args=false)
end
function convert_parameters(::Type{T}, p::MarkovPS{N}) where {T,N}
    return MarkovPS{N}(map(T, p.scales); check_args=false)
end

function check_frequency(p, ν)
    length(ν) == ndims(p) || throw(DimensionMismatch("Wrong frequency rank."))
    return nothing
end
# Specialize frequency tuple lengths in pointwise and broadcast evaluations.
function power_point(p::ConstantPowerSpectrum, ν::Vararg{Any,K}) where {K}
    check_frequency(p, ν)
    return p.power
end
function power_point(p::SinglePowerLaw, ν::Vararg{Any,K}) where {K}
    check_frequency(p, ν)
    r2, index, dc = promote(sum(abs2, ν), p.index, p.dc)

    # Arithmetic masks keep the DC expression finite during device differentiation.
    dc_mode = ifelse(iszero(r2), one(r2), zero(r2))
    r = sqrt(r2 + dc_mode)
    value = r^index
    return dc_mode * oftype(value, dc) + (one(dc_mode) - dc_mode) * value
end
function power_point(p::SaturatedPowerLaw, ν::Vararg{Any,K}) where {K}
    check_frequency(p, ν)
    r2 = sum(abs2, ν)
    r = sqrt(max(r2, p.inscale^2))
    return ifelse(r2 >= p.outscale^2, zero(r), r^p.index)
end
# Zero-power modes stay fixed at zero; use a zero derivative at this boundary.
function power_amplitude(power)
    zero_mode = ifelse(iszero(power), one(power), zero(power))
    return (one(power) - zero_mode) * sqrt(power + zero_mode)
end
function amplitude_point(p::AbstractPowerSpectrumModel, ν::Vararg{Any,K}) where {K}
    return power_amplitude(power_point(p, ν...))
end

struct SumPowerSpectrum{N,P<:Tuple} <: AbstractPowerSpectrumModel{N}
    models::P
end
function Base.:+(
    a::AbstractPowerSpectrumModel{N}, b::AbstractPowerSpectrumModel{N}
) where {N}
    models = promote_models((a, b))
    return SumPowerSpectrum{N,typeof(models)}(models)
end
power_point(p::SumPowerSpectrum, ν::Vararg{Any,K}) where {K} = summed_power(p.models, ν)
summed_power(models::Tuple{Any}, ν) = power_point(only(models), ν...)
function summed_power(models::Tuple, ν)
    return power_point(first(models), ν...) + summed_power(Base.tail(models), ν)
end

abstract type ModelModifier end
struct Renormalize{T} <: ModelModifier
    scale::T
    function Renormalize(s::T; check_args::Bool=true) where {T<:Number}
        (!check_args || (isfinite(s) && s >= 0)) ||
            throw(ArgumentError("Power scaling must be nonnegative."))
        return new{T}(s)
    end
end
struct Stretch{N,T} <: ModelModifier
    factors::NTuple{N,T}
end
function Stretch(s::Number...; check_args::Bool=true)
    f = promote(s...)
    !isempty(f) && (!check_args || all(v -> isfinite(v) && v > 0, f)) ||
        throw(ArgumentError("Stretch factors must be positive."))
    return Stretch(f)
end
struct Rotate{N,T} <: ModelModifier
    rows::NTuple{N,NTuple{N,T}}
end
function Rotate(angle::Number)
    angle isa Complex && throw(ArgumentError("Rotation angles must be real."))
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
function modify(p::AbstractPowerSpectrumModel{N}, ms::Vararg{ModelModifier,K}) where {N,K}
    for m in ms
        m isa Stretch &&
            length(m.factors) != N &&
            throw(DimensionMismatch("Stretch has the wrong rank."))
        m isa Rotate &&
            length(m.rows) != N &&
            throw(DimensionMismatch("Rotation has the wrong rank."))
    end
    return promoted_modified(
        p, ms, common_parameter_type((parameter_type(p), map(parameter_type, ms)...))
    )
end
function power_point(p::ModifiedPowerSpectrumModel, ν::Vararg{Any,K}) where {K}
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
    return modified_power(p, ms, map(row -> sum(row .* ν), m.rows))
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

# Promote built-in compositions before storing parameters. Mixed Float32/Float64
# fields can otherwise lose Enzyme adjoints on Julia 1.10 when sampled into Float32.
# Custom models retain their representation unless they opt into this internal trait.
parameter_type(::AbstractPowerSpectrumModel) = nothing
parameter_type(::ConstantPowerSpectrum{N,T}) where {N,T} = T
parameter_type(::SinglePowerLaw{N,T,S}) where {N,T,S} = promote_type(T, S)
parameter_type(::SaturatedPowerLaw{N,T}) where {N,T} = T
parameter_type(::Renormalize{T}) where {T} = T
parameter_type(::Stretch{N,T}) where {N,T} = T
parameter_type(::Rotate{N,T}) where {N,T} = T
function common_parameter_type(types)
    return any(isnothing, types) ? nothing : promote_type(types...)
end
parameter_type(p::SumPowerSpectrum) = common_parameter_type(map(parameter_type, p.models))
function parameter_type(p::ModifiedPowerSpectrumModel)
    return common_parameter_type((parameter_type(p.model), map(parameter_type, p.modifiers)...))
end
function convert_parameters(::Type{T}, p::ConstantPowerSpectrum{N}) where {T,N}
    return ConstantPowerSpectrum{N}(convert(T, p.power); check_args=false)
end
function convert_parameters(::Type{T}, p::SinglePowerLaw{N}) where {T,N}
    return SinglePowerLaw{N}(convert(T, p.index); dc=convert(T, p.dc), check_args=false)
end
function convert_parameters(::Type{T}, p::SaturatedPowerLaw{N}) where {T,N}
    return SaturatedPowerLaw{N,T}(
        convert(T, p.index), convert(T, p.inscale), convert(T, p.outscale)
    )
end
function convert_parameters(::Type{T}, m::Renormalize) where {T}
    return Renormalize(convert(T, m.scale); check_args=false)
end
function convert_parameters(::Type{T}, m::Stretch) where {T}
    return Stretch(map(v -> convert(T, v), m.factors))
end
function convert_parameters(::Type{T}, m::Rotate) where {T}
    return Rotate(map(row -> map(v -> convert(T, v), row), m.rows))
end
function convert_parameters(::Type{T}, p::SumPowerSpectrum{N}) where {T,N}
    models = map(m -> convert_parameters(T, m), p.models)
    return SumPowerSpectrum{N,typeof(models)}(models)
end
function raw_modified(p::AbstractPowerSpectrumModel{N}, ms) where {N}
    return ModifiedPowerSpectrumModel{N,typeof(p),typeof(ms)}(p, ms)
end
function convert_parameters(::Type{T}, p::ModifiedPowerSpectrumModel) where {T}
    return raw_modified(
        convert_parameters(T, p.model), map(m -> convert_parameters(T, m), p.modifiers)
    )
end
promoted_modified(p, ms, ::Nothing) = raw_modified(p, ms)
function promoted_modified(p, ms, ::Type{T}) where {T}
    return raw_modified(convert_parameters(T, p), map(m -> convert_parameters(T, m), ms))
end
function promote_models(models)
    T = common_parameter_type(map(parameter_type, models))
    return isnothing(T) ? models : map(m -> convert_parameters(T, m), models)
end

function frequency_axes(g::FieldGrid{N}, (::Type{T})=Float64) where {N,T}
    return ntuple(N) do d
        n = g.dims[d]
        [T(k <= (n - 1) ÷ 2 ? k : k - n) / T(n * g.spacing[d]) for k in 0:(n-1)]
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

@doc """
    power_point(model, frequencies...)

Evaluate nonnegative power at one frequency coordinate, in cycles per coordinate unit.
Custom spectrum types implement this method for their declared rank.
""" power_point

@doc """
    amplitude_point(model, frequencies...)

Evaluate the square root of pointwise power. Zero-power modes use a zero boundary
subgradient in transform execution; hard cutoffs remain nonsmooth.
""" amplitude_point

@doc """
    power_map(model, frequency_axes::Tuple)
    power_map(model, grid::FieldGrid)

Evaluate pointwise power on a Cartesian frequency grid. This does not perform the
conjugate-bin averaging used by SpectralScaling.
""" power_map

@doc """
    amplitude_map(model, grid)

Allocate the square root of a pointwise power map.
""" amplitude_map

@doc """
    power_map!(output, model, coordinates::Tuple)

Evaluate power using broadcastable coordinate arrays on the output backend. For a
2D grid, coordinate shapes are `(nx, 1)` and `(1, ny)`, not two unshaped axis vectors.
The caller owns and may reuse the coordinate arrays.
""" power_map!

@doc """
    amplitude_map!(output, model, coordinates::Tuple)

Write the square root of pointwise power into output using prepared coordinate arrays.
""" amplitude_map!

@doc """
    modify(model, modifiers...)

Compose spectral coordinate and power modifiers in the supplied order. The last
coordinate modifier acts first when evaluating the pullback to the base spectrum.
""" modify

@doc """
    renormed(model, factor)

Scale power by a finite nonnegative factor; signal amplitude changes by its square root.
""" renormed

@doc """
    stretched(model, factors...)

Evaluate the spectrum at frequencies divided by positive per-axis factors. Integrated
power is not renormalized.
""" stretched

@doc """
    rotated(model, angle_or_matrix)

Pull frequencies back through a 2D rotation angle or an N-dimensional proper orthogonal
matrix. Angles are in radians.
""" rotated

@doc """A nonnegative power multiplier. Disable checks only for validated dynamic parameterizations.""" Renormalize
@doc """Positive per-axis frequency widths. `check_args=false` supports compiled positive parameterizations.""" Stretch
@doc """A rotation of frequency coordinates, represented by fixed-size matrix rows.""" Rotate
