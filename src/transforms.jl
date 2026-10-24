"""A deterministic map with explicit input and output shapes."""
abstract type AbstractFieldTransform end
abstract type AbstractPreparedTransform end

"""A shape-preserving identity map."""
struct IdentityTransform{N} <: AbstractFieldTransform
    dims::NTuple{N, Int}
    function IdentityTransform(dims::NTuple{N, Integer}) where {N}
        return new{N}(size(FieldGrid(dims)))
    end
end
input_size(t::IdentityTransform) = t.dims
output_size(t::IdentityTransform) = t.dims
prepare(t::IdentityTransform, x) = (check_input(t, x); t)
apply(t::IdentityTransform, x) = (check_input(t, x); copy(x))

"""A true multidimensional, orthonormally normalized Hartley synthesis."""
struct HartleyTransform{N} <: AbstractFieldTransform
    dims::NTuple{N, Int}
    function HartleyTransform(dims::NTuple{N, Integer}) where {N}
        return new{N}(size(FieldGrid(dims)))
    end
end
input_size(t::HartleyTransform) = t.dims
output_size(t::HartleyTransform) = t.dims
struct PreparedHartley{N, P, T} <: AbstractPreparedTransform
    dims::NTuple{N, Int}
    plan::P
    scale::T
end
input_size(t::PreparedHartley) = t.dims
output_size(t::PreparedHartley) = t.dims
function prepare(t::HartleyTransform, x)
    check_input(t, x)
    p = plan_fht(x)
    return PreparedHartley(t.dims, p, inv(sqrt(eltype(x)(length(x)))))
end
function check_hartley_input(t::PreparedHartley, x)
    check_input(t, x)
    sample_eltype(x) == typeof(t.scale) ||
        throw(ArgumentError("Prepare a new Hartley plan for a different element type."))
    return nothing
end
apply(t::PreparedHartley, x) = (check_hartley_input(t, x); (t.plan * x) .* t.scale)
function apply!(y, t::PreparedHartley, x)
    check_buffers(y, t, x)
    check_hartley_input(t, x)
    Base.mightalias(y, x) &&
        throw(ArgumentError("Hartley input and output must not alias."))
    mul!(y, t.plan, x)
    y .*= t.scale
    return y
end

"""Multiply signal or coefficient amplitudes by a scalar or broadcastable array."""
struct ScaleTransform{N, S} <: AbstractFieldTransform
    dims::NTuple{N, Int}
    scale::S
    function ScaleTransform(dims::NTuple{N, Integer}, scale) where {N}
        shape = size(FieldGrid(dims))
        if scale isa AbstractArray
            ndims(scale) <= N && all(d -> size(scale, d) in (1, shape[d]), 1:N) ||
                throw(DimensionMismatch("Scale must broadcast to the transform shape."))
        end
        return new{N, typeof(scale)}(shape, scale)
    end
end
input_size(t::ScaleTransform) = t.dims
output_size(t::ScaleTransform) = t.dims
prepare(t::ScaleTransform, x) = (check_input(t, x); t)
apply(t::ScaleTransform, x) = (check_input(t, x); sample_eltype(x).(x .* t.scale))

"""
    SpectralScaling(grid, model)

Multiply Hartley coefficients by the square root of discrete power.
Conjugate frequency pairs are averaged on the sampled grid, including Nyquist aliases.
The power model is reevaluated on every application, so parameter changes are visible.
"""
struct SpectralScaling{G, P} <: AbstractFieldTransform
    grid::G
    model::P
    function SpectralScaling(grid::FieldGrid, model::AbstractPowerSpectrumModel)
        ndims(grid) == ndims(model) || throw(DimensionMismatch("Spectrum rank mismatch."))
        return new{typeof(grid), typeof(model)}(grid, model)
    end
end
input_size(t::SpectralScaling) = size(t.grid)
output_size(t::SpectralScaling) = size(t.grid)
struct PreparedSpectralScaling{N, P, C, R} <: AbstractPreparedTransform
    dims::NTuple{N, Int}
    model::P
    coordinates::C
    reflected::R
end
input_size(t::PreparedSpectralScaling) = t.dims
output_size(t::PreparedSpectralScaling) = t.dims
function prepare(t::SpectralScaling, x)
    check_input(t, x)
    coords = frequency_coordinates(t.grid, x)
    reflected = frequency_coordinates(t.grid, x; reflected = true)
    return PreparedSpectralScaling(input_size(t), t.model, coords, reflected)
end
function paired_power(model, ν::Vararg{Any, K}) where {K}
    n = ndims(model)
    p = power_point(model, ν[1:n]...)
    q = power_point(model, ν[(n + 1):(2n)]...)
    return (p + q) / 2
end
function paired_amplitude(model, ν::Vararg{Any, K}) where {K}
    return power_amplitude(paired_power(model, ν...))
end
function spectral_amplitude(t::PreparedSpectralScaling, model)
    return paired_amplitude.(Ref(model), t.coordinates..., t.reflected...)
end
apply(t::PreparedSpectralScaling, x) = apply(t, x, t.model)
function apply(t::PreparedSpectralScaling, x, model::AbstractPowerSpectrumModel)
    check_input(t, x)
    ndims(model) == length(t.dims) || throw(DimensionMismatch("Spectrum rank mismatch."))
    return sample_eltype(x).(x .* spectral_amplitude(t, model))
end
function apply!(y, t::PreparedSpectralScaling, x, model::AbstractPowerSpectrumModel = t.model)
    check_buffers(y, t, x)
    y .= x .* spectral_amplitude(t, model)
    return y
end

"""A left-to-right composition: `TransformChain(a, b)` computes `b(a(x))`."""
struct TransformChain{T <: Tuple} <: AbstractFieldTransform
    transforms::T
end
function TransformChain(ts::Union{AbstractFieldTransform, AbstractPreparedTransform}...)
    isempty(ts) && throw(ArgumentError("A transform chain cannot be empty."))
    for i in 2:length(ts)
        output_size(ts[i - 1]) == input_size(ts[i]) ||
            throw(DimensionMismatch("Transform shapes do not compose."))
    end
    return TransformChain(ts)
end
input_size(t::TransformChain) = input_size(first(t.transforms))
output_size(t::TransformChain) = output_size(last(t.transforms))
function prepare(t::TransformChain, x)
    check_input(t, x)
    return TransformChain(
        map(t.transforms) do part
            prototype = similar(x, eltype(x), input_size(part))
            prepare(part, prototype)
        end,
    )
end
apply(t::TransformChain, x) = apply_chain(t.transforms, x)
apply_chain(::Tuple{}, x) = x
apply_chain(ts::Tuple, x) = apply_chain(Base.tail(ts), apply(first(ts), x))

prepare(t::AbstractPreparedTransform, x) = (check_input(t, x); t)
apply(t::AbstractFieldTransform, x) = apply(prepare(t, x), x)
function apply!(y, t::Union{AbstractFieldTransform, AbstractPreparedTransform}, x)
    check_buffers(y, t, x)
    copyto!(y, apply(t, x))
    return y
end
function check_input(t, x)
    size(x) == input_size(t) || throw(DimensionMismatch("Wrong transform input shape."))
    eltype(x) <: AbstractFloat ||
        throw(ArgumentError("Transforms require floating-point arrays."))
    return nothing
end
function check_buffers(y, t, x)
    check_input(t, x)
    size(y) == output_size(t) || throw(DimensionMismatch("Wrong transform output shape."))
    eltype(y) == eltype(x) ||
        throw(ArgumentError("Input and output element types must match."))
    return nothing
end

sample_eltype(x::AbstractArray) = eltype(x)

@doc """
    input_size(transform)

Return the shape of one latent/coefficient input, excluding batch axes.
""" input_size

@doc """
    output_size(transform)

Return the shape of one transformed signal, excluding batch axes.
""" output_size

@doc """
    prepare(transform, prototype)
    prepare(field, latent_prototype)

Prepare fixed-shape execution metadata for the prototype's backend and element type.
Call outside differentiation and compilation. Prepare each original descriptor separately
for concurrent tasks. Model-dependent spectral amplitudes are not cached.
""" prepare

@doc """
    apply(transform, x)
    apply(prepared_spectral, x, model)

Apply a deterministic transform. An unprepared descriptor is prepared automatically;
reuse a prepared object for repeated execution. Built-in transforms preserve the input
floating-point element type. The three-argument form reevaluates a new spectrum model.
""" apply

@doc """
    apply!(output, transform, input)

Write a transformed signal to a matching preallocated output and return it. Internal
allocations depend on the backend; this is not an allocation-free guarantee. Hartley
execution rejects overlapping input and output storage.
""" apply!
