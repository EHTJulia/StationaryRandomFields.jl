export stretched, rotated, renormed, modify, Stretch, Renormalize, Rotate

"""
    $(TYPEDEF)

General type for a model modifier. These transform any model
using simple Fourier transform properties. To modify a model
you can use the [`ModifiedPowerSpectrumModel{N}`](@ref) constructor or the [`modify`](@ref)
function.

```julia-repl
julia> visanalytic(stretched(Disk(), 2.0, 2.0)) == visanalytic(Disk())
true
```

To implement a model transform you need to specify the following methods:
- [`transform_ν`](@ref)
- [`transform_image`](@ref)
- [`scale_fourier`](@ref)
- [`scale_image`](@ref)
- [`radialextent`](@ref)
See thee docstrings of those methods for guidance on implementation details.

Additionally these methods assume the modifiers are of the form

I(x,y) -> fᵢ(x,y)I(gᵢ(x,y))
V(u,v) -> fᵥ(u,v)V(gᵥ(u,v))

where `g` are the transform_image/uv functions and `f` are the scale_image/uv
function.

"""
abstract type ModelModifier{T} end

"""
    scale_fourier(model::AbstractModifier, ν...)

Returns a number on how to scale the fourier coefficients at the frequency coordinate `ν` for an modified `model`
"""
function scale_fourier end


"""
    transform_ν(model::AbstractModifier, ν...)

Returns a transformed frequency coordinate `ν` according to the `model` modifier
"""
function transform_ν end

unitscale(T) = one(T)
unitscale(T, M) = one(T)

"""
    $(TYPEDEF)

Container type for models that have been transformed in some way.
For a list of potential modifiers or transforms see `subtypes(ModelModifiers)`.

# Fields
$(FIELDS)
"""
struct ModifiedPowerSpectrumModel{N,M<:AbstractPowerSpectrumModel{N},T<:Tuple} <: AbstractPowerSpectrumModel{N}
    """base model"""
    model::M
    """model transforms"""
    transform::T
end

function Base.show(io::IO, mi::ModifiedPowerSpectrumModel)
    #io = IOContext(io, :compact=>true)
    #s = summary(mi)
    println(io, "ModifiedPowerSpectrumModel")
    println(io, "  base model: ", summary(mi.model))
    println(io, "  Modifiers:")
    for i in eachindex(mi.transform)
        println(io, "    $i. ", summary(mi.transform[i]))
    end
end


"""
    unmodified(model::ModifiedPowerSpectrumModel)

Returns the un-modified model

### Example
```julia-repl
julia> m = stretched(rotated(Gaussian(), π/4), 2.0, 1.0)
julia> umodified(m) == Gaussian()
true
```
"""
unmodified(model::ModifiedPowerSpectrumModel) = model.model
unmodified(model::AbstractPowerSpectrumModel) = model

"""
    basemodel(model::ModifiedPowerSpectrumModel)

Returns the ModifiedPowerSpectrumModel with the last transformation stripped.

# Example
```julia-repl
julia> basemodel(stretched(Disk(), 1.0, 2.0)) == Disk()
true
```
"""
basemodel(model::ModifiedPowerSpectrumModel) = ModifiedPowerSpectrumModel(model.model, Base.front(model.transform))
basemodel(model::ModifiedPowerSpectrumModel{N,M,Tuple{}}) where {N,M} = model

@inline fourieranalytic(::Type{ModifiedPowerSpectrumModel{N,M,T}}) where {N,M,T} = fourieranalytic(M)

@inline function ModifiedPowerSpectrumModel(m::AbstractPowerSpectrumModel, t::ModelModifier)
    return ModifiedPowerSpectrumModel(m, (t,))
end

@inline function ModifiedPowerSpectrumModel(m::ModifiedPowerSpectrumModel, t::ModelModifier)
    model = m.model
    t0 = m.transform
    return ModifiedPowerSpectrumModel(model, (t0..., t))
end

function modify_fourier(model, t::Tuple, scale, ν...)
    νt = transform_ν(model, last(t), ν...)
    scalet = scale_fourier(model, last(t), ν...)
    return modify_fourier(model, Base.front(t), scalet * scale, νt...)
end
modify_fourier(model, ::Tuple{}, scale, ν...) = scale, ν

"""
    modify(m::AbstractPowerSpectrumModel, transforms...)

Modify a given `model` using the set of `transforms`. This is the most general
function that allows you to apply a sequence of model transformation for example

```julia-repl
modify(Gaussian(), Stretch(2.0, 1.0), Rotate(π/4), Shift(1.0, 2.0), Renorm(2.0))
```
will create a asymmetric Gaussian with position angle `π/4` shifted to the position
(1.0, 2.0) with a flux of 2 Jy. This is similar to Flux's chain.
"""
function modify(m::AbstractPowerSpectrumModel, transforms...)
    return ModifiedPowerSpectrumModel(m, transforms)
end

@inline function fourier_point(m::ModifiedPowerSpectrumModel, ν...)
    mbase = m.model
    transform = m.transform
    scale, ν = modify_fourier(mbase, transform, unitscale(Complex{eltype(u)}, typeof(mbase)), ν...)
    scale * fourier_point(mbase, ν...)
end

function modify_fourier(model, transform::Tuple, scale, ν::AbstractVector...)
    uvscale = modify_fourier.(Ref(model), Ref(transform), scale, ν::AbstractVector...)
    return first.(uvscale), last.(uvscale)
end

function __extract_tangent(dm::ModifiedPowerSpectrumModel)
    tm = __extract_tangent(dm.model)
    dtm = dm.transform
    ttm = map(x -> Tangent{typeof(x)}(; ntfromstruct(x)...), dtm)
    tm = Tangent{typeof(dm)}(model=tm, transform=ttm)
end

"""
    $(TYPEDEF)

Renormalizes the flux of the model to the new value `scale*flux(model)`.
We have also overloaded the Base.:* operator as syntactic sugar
although I may get rid of this.


# Example
```julia-repl
julia> modify(Gaussian(), Renormalize(2.0)) == 2.0*Gaussian()
true
```
"""
struct Renormalize{T} <: ModelModifier{T}
    scale::T
end


"""
    $(SIGNATURES)

Renormalizes the model `m` to have total flux `f*flux(m)`.
This can also be done directly by calling `Base.:*` i.e.,

```julia-repl
julia> renormed(m, f) == f*M
true
```
"""
renormed(model::M, f) where {M<:AbstractPowerSpectrumModel} = ModifiedPowerSpectrumModel(model, Renormalize(f))
Base.:*(model::AbstractPowerSpectrumModel, f::Number) = renormed(model, f)
Base.:*(f::Number, model::AbstractPowerSpectrumModel) = renormed(model, f)
Base.:/(f::Number, model::AbstractPowerSpectrumModel) = renormed(model, inv(f))
Base.:/(model::AbstractPowerSpectrumModel, f::Number) = renormed(model, inv(f))
# Dispatch on RenormalizedModel so that I just make a new RenormalizedModel with a different f
# This will make it easier on the compiler.
# Base.:*(model::ModifiedPowerSpectrumModel{N}, f::Number) = renormed(model.model, model.scale*f)
# Overload the unary negation operator to be the same model with negative flux
Base.:-(model::AbstractPowerSpectrumModel) = renormed(model, -1)

@inline transform_ν(m, ::Renormalize, ν...) = (ν...,)
@inline scale_fourier(::M, transform::Renormalize, ν...) where {M} = transform.scale * unitscale(typeof(transform.scale), M)

"""
    Stretch(α, β)
    Stretch(r)

Stretched the model in the x and y directions, i.e. the new intensity is
    Iₛ(x,y) = 1/(αβ) I(x/α, y/β),
where were renormalize the intensity to preserve the models flux.

# Example
```julia-repl
julia> modify(Gaussian(), Stretch(2.0)) == stretched(Gaussian(), 2.0, 1.0)
true
```

If only a single argument is given it assumes the same stretch is applied in both direction.

```julia-repl
julia> Stretch(2.0) == Stretch(2.0, 2.0)
true
```
"""
struct Stretch{T,N} <: ModelModifier{T}
    α::NTuple{N,T}
    function Stretch(scale::T...) where {T}
        new{T,length(scale)}(scale)
    end
end

"""
    $(SIGNATURES)

Stretches the model `m` according to the formula
    Iₛ(x,y) = 1/(αβ) I(x/α, y/β),
where were renormalize the intensity to preserve the models flux.
"""
stretched(model, α...) = ModifiedPowerSpectrumModel{N}(model, Stretch(α...))

@inline transform_ν(m, transform::Stretch, ν...) = ν .* transform.α
@inline scale_fourier(::M, ::Stretch{T}, ν...) where {M,T} = unitscale(T, M)


"""
    Rotate(ξ)

Type for the rotated model. This is more fine grained constrol of
rotated model.

# Example
```julia-repl
julia> modify(Gaussian(), Rotate(2.0)) == rotated(Gaussian(), 2.0)
true
```
"""
struct Rotate{T} <: ModelModifier{T}
    s::T
    c::T
    function Rotate(ξ::F) where {F}
        s, c = sincos(ξ)
        return new{F}(s, c)
    end
end


"""
    $(SIGNATURES)

Rotates the model by an amount `ξ` in radians in the clockwise direction.
"""
rotated(model, ξ) = ModifiedPowerSpectrumModel(model, Rotate(ξ))

"""
    $(SIGNATURES)

Returns the rotation angle of the rotated `model`
"""
posangle(model::Rotate) = atan(model.s, model.c)

@inline function transform_ν(m, transform::Rotate, ν...)
    if length(ν) == 1
        return ν
    else
        s, c = transform.s, transform.c
        return (c * u - s * v, +s * u + c * v, ν[3:end])
    end
end

@inline scale_fourier(::Rotate{T}, ν...) where {T} = one(T)
