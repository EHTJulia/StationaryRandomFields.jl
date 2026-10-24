module StationaryRandomFieldsCUDAExt

using CUDA
using Random
using Distributions
using StationaryRandomFields
using Wavelets
import Distributions: _rand!
import StationaryRandomFields: apply, wavelet_analysis

function apply(t::StationaryRandomFields.PreparedWavelet, x::CuArray)
    StationaryRandomFields.check_input(t, x)
    y = similar(x)
    if ndims(x) <= 3 && !isnothing(Base.get_extension(Wavelets, :WaveletsGPUExt))
        Wavelets.idwt!(y, x, t.wavelet, t.levels)
    else
        copyto!(y, StationaryRandomFields.wavelet_synthesis(x, t.filters, t.levels))
    end
    return y
end

function wavelet_analysis(t::StationaryRandomFields.PreparedWavelet, x::CuArray)
    StationaryRandomFields.check_input(t, x)
    y = similar(x)
    if ndims(x) <= 3 && !isnothing(Base.get_extension(Wavelets, :WaveletsGPUExt))
        Wavelets.dwt!(y, x, t.wavelet, t.levels)
    else
        copyto!(y, StationaryRandomFields.wavelet_analysis(x, t.filters, t.levels))
    end
    return y
end

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
