"""
    TransformedRandomField(latent, transform)

A signal sampler defined by `x = apply(transform, z)` with `z ~ latent`.
Only the latent distribution provides `logpdf`; signal density is deliberately unavailable.
"""
struct TransformedRandomField{N,L,T} <: AbstractStationaryRandomField{N}
    latent::L
    transform::T
end
function TransformedRandomField(latent::AbstractStationaryRandomField, transform)
    size(latent) == input_size(transform) ||
        throw(DimensionMismatch("Latent shape mismatch."))
    n = length(output_size(transform))
    return TransformedRandomField{n,typeof(latent),typeof(transform)}(latent, transform)
end
Base.size(f::TransformedRandomField) = output_size(f.transform)
Base.eltype(f::TransformedRandomField) = eltype(f.latent)
Base.eltype(::Type{TransformedRandomField{N,L,T}}) where {N,L,T} = eltype(L)

function signal_density_error()
    throw(
        ArgumentError(
            "Signal logpdf is not provided. Evaluate logpdf(field.latent, z) in latent space.",
        ),
    )
end
function Distributions._logpdf(::TransformedRandomField, ::AbstractArray{<:Real})
    signal_density_error()
end
function Distributions._logpdf(::TransformedRandomField{2}, ::AbstractMatrix{<:Real})
    signal_density_error()
end
function Distributions.insupport(::TransformedRandomField, ::Any)
    throw(
        ArgumentError(
            "Signal support membership is not implemented for general transforms."
        ),
    )
end

function prepare(f::TransformedRandomField, prototype)
    size(prototype) == size(f.latent) ||
        throw(DimensionMismatch("Latent prototype mismatch."))
    return TransformedRandomField(f.latent, prepare(f.transform, prototype))
end
function Distributions._rand!(
    rng::AbstractRNG, f::TransformedRandomField{N}, x::AbstractArray{<:Real,N}
) where {N}
    z = similar(x, eltype(x), size(f.latent))
    rand!(rng, f.latent, z)
    return apply!(x, f.transform, z)
end
