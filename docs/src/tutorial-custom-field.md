# Custom fields

For independent coefficients, pass any continuous scalar distribution to
`UnivariateRandomField`. This example uses a Student-t prior:

```@example custom_field
using StationaryRandomFields, Distributions, Random

prior = UnivariateRandomField((32, 32), TDist(5.0))
z = rand(MersenneTwister(12), prior)
logpdf(prior, z)
```

For a new field type, implement its shape, element type, sampling, support, and
log density. Here is a one-dimensional field of independent uniform coefficients:

```@example custom_field
struct UnitUniformField <: AbstractStationaryRandomField{1}
    n::Int
    function UnitUniformField(n::Integer)
        n > 0 || throw(ArgumentError("Length must be positive."))
        return new(Int(n))
    end
end

Base.size(f::UnitUniformField) = (f.n,)
Base.eltype(::Type{UnitUniformField}) = Float64
Base.eltype(::UnitUniformField) = Float64
Distributions.insupport(f::UnitUniformField, z::AbstractVector{<:Real}) =
    size(z) == size(f) && all(x -> 0 <= x <= 1, z)
Distributions._logpdf(f::UnitUniformField, z::AbstractVector{<:Real}) =
    insupport(f, z) ? 0.0 : -Inf
function Distributions._rand!(rng::AbstractRNG, f::UnitUniformField, z::AbstractVector{<:Real})
    size(z) == size(f) || throw(DimensionMismatch("Wrong field shape."))
    return rand!(rng, z)
end

prior = UnitUniformField(32)
z = rand(MersenneTwister(12), prior)
(logpdf(prior, z), size(z))
```

The field can then be combined with a transform:

```@example custom_field
transform = prepare(HartleyTransform(), prior)
field = TransformedRandomField(prior, transform)
x = rand(MersenneTwister(13), field)
size(x)
```

```@example custom_field
@assert logpdf(prior, z) == 0.0 # hide
@assert logpdf(prior, fill(-1.0, 32)) == -Inf # hide
@assert x ≈ apply(transform, rand(MersenneTwister(13), prior)) # hide
nothing # hide
```

This sampling method uses a CPU RNG. Custom fields need backend-specific sampling
methods to generate coefficients on a GPU. Prior evaluation belongs to the
[latent field](math-latent.md).
