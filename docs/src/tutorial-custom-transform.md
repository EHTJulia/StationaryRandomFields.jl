# Custom transforms

A transform descriptor defines the operation; preparation resolves its shape.
This example adds a constant offset to a signal:

```@example custom_transform
using StationaryRandomFields, Random

struct OffsetTransform{T} <: AbstractFieldTransform
    offset::T
end

struct PreparedOffset{N, T} <: StationaryRandomFields.AbstractPreparedTransform
    dims::NTuple{N, Int}
    offset::T
end
StationaryRandomFields.input_size(t::PreparedOffset) = t.dims
StationaryRandomFields.output_size(t::PreparedOffset) = t.dims
function StationaryRandomFields.prepare(t::OffsetTransform, x::AbstractArray)
    return PreparedOffset(size(x), t.offset)
end
function StationaryRandomFields.apply(t::PreparedOffset, x)
    size(x) == t.dims || throw(DimensionMismatch("Wrong input shape."))
    return x .+ t.offset
end
```

The descriptor stores only the offset. `prepare` resolves the input shape and returns
an execution object. Define `input_size`, `output_size`, and `apply` on that object.
The array-specific preparation method also works through `prepare(transform, prior)`.

```@example custom_transform
prior = UnivariateGaussianRandomField((32, 32))
z = rand(MersenneTwister(12), prior)
transform = TransformChain(
    HartleyTransform(), OffsetTransform(2.0)
)
prepared = prepare(transform, prior)
x = apply(prepared, z)
field = TransformedRandomField(prior, prepared)
size(rand(field))
```

```@example custom_transform
@assert x ≈ apply(prepare(HartleyTransform(), prior), z) .+ 2.0 # hide
@assert rand(MersenneTwister(12), field) ≈ x # hide
nothing # hide
```

`TransformChain` applies its members in order. A specialized `apply!` is optional;
the default copies the result of `apply` into the supplied output buffer.
