module StationaryRandomFieldsReactantExt

using Reactant
using Random
using Distributions
using StationaryRandomFields
import Distributions: _rand!, logpdf
import StationaryRandomFields: check_input, latent_logpdf, input_size, sample_eltype

const ReactantArray{T, N} = Union{
    Reactant.AnyTracedRArray{T, N}, Reactant.AnyConcreteRArray{T, N},
}

# Traced scalar element types do not subtype Real, so Distributions' host
# array fallbacks do not dispatch on traced arrays.
function Distributions.logpdf(
        f::UnivariateRandomField{N}, x::ReactantArray{T, N}
    ) where {T, N}
    size(x) == size(f) || throw(DimensionMismatch("Wrong latent shape."))
    return latent_logpdf(f.dist, x)
end
function Distributions.logpdf(
        f::UnivariateRandomField{N, D},
        x::ReactantArray{T, N};
        location = mean(f.dist),
        scale = StationaryRandomFields.prior_scale(f.dist),
        check_args = true,
    ) where {T, N, D <: Union{Normal, Laplace}}
    return StationaryRandomFields.parameterized_logpdf(f, x, location, scale, check_args)
end
sample_eltype(::ReactantArray{T}) where {T} = T

function check_input(t, x::ReactantArray{T}) where {T}
    size(x) == input_size(t) || throw(DimensionMismatch("Wrong transform input shape."))
    T <: AbstractFloat || throw(ArgumentError("Transforms require floating-point arrays."))
    return nothing
end

function _rand!(
        rng::AbstractRNG, f::UnivariateRandomField{N, D}, x::ReactantArray{T, N}
    ) where {N, D <: Normal, T <: AbstractFloat}
    Random.randn!(rng, x)
    x .= T(mean(f.dist)) .+ T(std(f.dist)) .* x
    return x
end
function _rand!(
        rng::AbstractRNG, f::UnivariateRandomField{N, D}, x::ReactantArray{T, N}
    ) where {N, D <: Laplace, T <: AbstractFloat}
    Random.rand!(rng, x)
    x .= clamp.(x, eps(T), one(T) - eps(T)) .- T(0.5)
    x .= T(f.dist.μ) .- T(f.dist.θ) .* sign.(x) .* log1p.(-2 .* abs.(x))
    return x
end

function _rand!(
        rng::AbstractRNG, f::TransformedRandomField{N}, x::ReactantArray{T, N}
    ) where {T, N}
    size(x) == size(f) || throw(DimensionMismatch("Wrong signal shape."))
    z = similar(x, T, size(f.latent))
    rand!(rng, f.latent, z)
    return apply!(x, f.transform, z)
end

end
