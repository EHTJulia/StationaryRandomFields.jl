"""A transform description whose grid is resolved by `prepare`."""
struct TransformSpec{K, P} <: AbstractFieldTransform
    parameters::P
end
ScaleTransform(scale) = TransformSpec{:scale}(scale)
HartleyTransform() = TransformSpec{:hartley}(nothing)
TransformSpec{K}(p) where {K} = TransformSpec{K, typeof(p)}(p)
"""
    SpectralScaling(model; normalization = :none)

Describe square-root power scaling without fixing a grid. `prepare` resolves the
coordinates. `:unit_variance` divides effective power by its grid mean on every
application, including after model updates. It assumes unit-variance latent
coefficients and requires positive total power; it does not normalize each draw.
"""
function SpectralScaling(model::AbstractPowerSpectrumModel; normalization = :none)
    normalization in (:none, :unit_variance) || throw(ArgumentError("Unknown spectral normalization."))
    return TransformSpec{:spectral}((; model, normalization))
end
function WaveletTransform(wt = Wavelets.wavelet(Wavelets.WT.haar); levels = nothing)
    return TransformSpec{:wavelet}((; wavelet = wt, levels))
end
has_shape(t::AbstractFieldTransform) = applicable(input_size, t)
has_shape(::AbstractPreparedTransform) = true
has_shape(::TransformSpec) = false
has_shape(t::TransformChain) = all(has_shape, t.transforms)

"""Resolve shape and coordinates from a latent field before preparing execution metadata.

`prototype` selects array storage without using its values. The default is a CPU
array with the latent field's shape and element type. Specify either `grid` or
`spacing`; omitted spacing is one along every axis. Density remains `logpdf(latent, z)`.
"""
function prepare(
        t::AbstractFieldTransform, latent::AbstractStationaryRandomField;
        prototype = nothing, grid = nothing, spacing = nothing
    )
    x = isnothing(prototype) ? Array{eltype(latent)}(undef, size(latent)) : prototype
    size(x) == size(latent) || throw(DimensionMismatch("Prototype and latent shapes differ."))
    sample_eltype(x) == eltype(latent) || throw(ArgumentError("Prototype and latent element types differ."))
    return prepare_on_grid(t, x, resolve_grid(x, grid, spacing))
end
function resolve_grid(x, grid, spacing)
    isnothing(grid) || isnothing(spacing) || throw(ArgumentError("Specify grid or spacing, not both."))
    g = isnothing(grid) ? (isnothing(spacing) ? FieldGrid(size(x)) : FieldGrid(size(x); spacing)) : grid
    size(g) == size(x) || throw(DimensionMismatch("Grid and input shapes differ."))
    return g
end
prepare_on_grid(t, x, grid) = prepare(t, x)
prepare_on_grid(t::TransformSpec{:scale}, x, grid) = prepare(ScaleTransform(size(grid), t.parameters), x)
prepare_on_grid(::TransformSpec{:hartley}, x, grid) = prepare(HartleyTransform(size(grid)), x)
function prepare_on_grid(t::TransformSpec{:wavelet}, x, grid)
    return prepare(WaveletTransform(size(grid), t.parameters.wavelet; levels = t.parameters.levels), x)
end
struct UnitVarianceScaling{P} <: AbstractPreparedTransform
    spectral::P
end
input_size(t::UnitVarianceScaling) = input_size(t.spectral)
output_size(t::UnitVarianceScaling) = output_size(t.spectral)
function apply(t::UnitVarianceScaling, x, model = t.spectral.model)
    check_input(t, x)
    a = spectral_amplitude(t.spectral, model)
    return sample_eltype(x).(x .* a ./ sqrt(sum(abs2, a) / length(a)))
end
function prepare_on_grid(t::TransformSpec{:spectral}, x, grid)
    p = prepare(SpectralScaling(grid, t.parameters.model), x)
    return t.parameters.normalization === :none ? p : UnitVarianceScaling(p)
end
function prepare(t::TransformSpec, x::AbstractArray; grid = nothing, spacing = nothing)
    return prepare_on_grid(t, x, resolve_grid(x, grid, spacing))
end
function prepare_on_grid(t::TransformChain, x, grid)
    current = x
    parts = map(t.transforms) do part
        prepared = prepare_on_grid(part, current, grid)
        current = similar(x, sample_eltype(x), output_size(prepared))
        prepared
    end
    return TransformChain(parts...)
end
# Resolve this intersection explicitly to keep latent preparation unambiguous.
function prepare(t::TransformChain, latent::AbstractStationaryRandomField; kwargs...)
    return invoke(prepare, Tuple{AbstractFieldTransform, AbstractStationaryRandomField}, t, latent; kwargs...)
end
function Base.show(io::IO, t::TransformSpec{K}) where {K}
    return if K === :scale
        print(io, "ScaleTransform(scale=")
        t.parameters isa AbstractArray ? summary(io, t.parameters) : show(io, t.parameters)
        print(io, ")")
    elseif K === :hartley
        print(io, "HartleyTransform(normalization=orthonormal)")
    elseif K === :spectral
        print(io, "SpectralScaling(model=")
        show(io, t.parameters.model)
        print(io, ", action=√P scaling, normalization=", t.parameters.normalization, ")")
    else
        print(
            io, "WaveletTransform(basis=", Wavelets.WT.name(t.parameters.wavelet),
            ", levels=", isnothing(t.parameters.levels) ? "maximum" : t.parameters.levels,
            ", boundary=periodic, direction=synthesis)"
        )
    end
end
function Base.show(io::IO, t::UnitVarianceScaling)
    print(io, "UnitVarianceScaling(")
    show(io, t.spectral)
    return print(io, ")")
end

for T in (IdentityTransform, HartleyTransform, ScaleTransform, SpectralScaling, WaveletTransform)
    @eval function prepare(t::$T, latent::AbstractStationaryRandomField; kwargs...)
        return invoke(prepare, Tuple{AbstractFieldTransform, AbstractStationaryRandomField}, t, latent; kwargs...)
    end
end
