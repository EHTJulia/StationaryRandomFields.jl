module StationaryRandomFieldsCUDAExt

using CUDA
using Random
using Distributions
using StationaryRandomFields
import Distributions: _rand!
import StationaryRandomFields: _rand_normal_field!, _rand_laplace_field!

function _rand!(
        rng::AbstractRNG, f::UnivariateRandomField{N, D}, x::CuArray{T, N}
    ) where {N, D <: Normal, T <: AbstractFloat}
    return _rand_normal_field!(rng, f, x)
end

function _rand!(
        rng::AbstractRNG, f::UnivariateRandomField{N, D}, x::CuArray{T, N}
    ) where {N, D <: Laplace, T <: AbstractFloat}
    return _rand_laplace_field!(rng, f, x)
end

Random.rand!(f::AbstractStationaryRandomField, x::CuArray) = rand!(CUDA.default_rng(), f, x)

end
