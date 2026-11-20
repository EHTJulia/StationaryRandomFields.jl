module StationaryRandomFieldsCUDAExt

using CUDA
using Random
using Distributions
using StationaryRandomFields
using Wavelets
import Distributions: _rand!
import StationaryRandomFields: _rand_normal_field!, _rand_laplace_field!, apply, wavelet_analysis

# The shared synthesis previously inferred `Any` for CUDA arrays in 1D and 2D.
# Keep the native Wavelets.jl GPU path for its supported ranks and retain the
# shared implementation as the fallback for higher-rank arrays.
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
        rng::AbstractRNG, f::IndependentRandomField{N, D}, x::CuArray{T, N}
    ) where {N, D <: Normal, T <: AbstractFloat}
    return _rand_normal_field!(rng, f, x)
end

function _rand!(
        rng::AbstractRNG, f::IndependentRandomField{N, D}, x::CuArray{T, N}
    ) where {N, D <: Laplace, T <: AbstractFloat}
    return _rand_laplace_field!(rng, f, x)
end

Random.rand!(f::AbstractStationaryRandomField, x::CuArray) = rand!(CUDA.default_rng(), f, x)

end
