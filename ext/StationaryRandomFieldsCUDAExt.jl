module StationaryRandomFieldsCUDAExt

using CUDA
using Random
using Distributions
using StationaryRandomFields
import Distributions: _rand!

function _rand!(
        rng::AbstractRNG, f::UnivariateRandomField{N, D}, x::CuArray{T, N}
    ) where {N, D <: Normal, T <: AbstractFloat}
    Random.randn!(rng, x)
    x .= T(mean(f.dist)) .+ T(std(f.dist)) .* x
    return x
end
function _rand!(
        rng::AbstractRNG, f::UnivariateRandomField{N, D}, x::CuArray{T, N}
    ) where {N, D <: Laplace, T <: AbstractFloat}
    Random.rand!(rng, x)
    # Map the midpoint of each representable uniform bin into (0, 1).
    x .= clamp.(x, eps(T), one(T) - eps(T)) .- T(0.5)
    x .= T(f.dist.μ) .- T(f.dist.θ) .* sign.(x) .* log1p.(-2 .* abs.(x))
    return x
end
Random.rand!(f::AbstractStationaryRandomField, x::CuArray) = rand!(CUDA.default_rng(), f, x)

end
