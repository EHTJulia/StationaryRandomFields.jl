"""
    WaveletTransform(dims, wavelet = wavelet(WT.haar); levels = nothing)

Periodic, orthogonal, multilevel wavelet synthesis with Wavelets.jl-compatible
coefficient packing.

The transform uses a separable periodic filter bank along each array axis,
supporting arbitrary dimensions on CPU, CUDA, and Reactant backends.
Wavelets.jl provides wavelet definitions and filter coefficients.

By default, use the largest level count for which every dimension is divisible
by `2^levels`. Pass an explicit integer to use fewer levels. An odd-sized axis
limits the default to zero levels (the identity).

Each dimension must be divisible by `2^levels`.
"""
struct WaveletTransform{N, W} <: AbstractFieldTransform
    dims::NTuple{N, Int}
    wavelet::W
    levels::Int
end
function WaveletTransform(
        dims::NTuple{N, Integer}, wt = Wavelets.wavelet(Wavelets.WT.haar); levels::Union{Nothing, Integer} = nothing
    ) where {N}
    wt isa Wavelets.WT.OrthoFilter{Wavelets.WT.PerBoundary} ||
        throw(ArgumentError("Only periodic orthogonal filter wavelets are supported."))
    shape = size(FieldGrid(dims))
    levels = isnothing(levels) ? minimum(trailing_zeros, shape) : levels
    levels >= 0 && levels <= minimum(trailing_zeros, shape) ||
        throw(ArgumentError("Each dimension must be divisible by 2^levels."))
    return WaveletTransform(shape, wt, Int(levels))
end
input_size(t::WaveletTransform) = t.dims
output_size(t::WaveletTransform) = t.dims
struct PreparedWavelet{N, W, F} <: AbstractPreparedTransform
    dims::NTuple{N, Int}
    wavelet::W
    filters::F
    levels::Int
end
input_size(t::PreparedWavelet) = t.dims
output_size(t::PreparedWavelet) = t.dims
function prepare(t::WaveletTransform, x)
    check_input(t, x)
    return PreparedWavelet(
        t.dims, deepcopy(t.wavelet), Tuple(eltype(x).(t.wavelet.qmf)), t.levels
    )
end
function apply(t::PreparedWavelet, x::AbstractArray)
    check_input(t, x)
    return wavelet_synthesis(x, t.filters, t.levels)
end
function wavelet_analysis(t::PreparedWavelet, x::AbstractArray)
    return wavelet_analysis(x, t.filters, t.levels)
end

function axis_slice(x, d, range)
    indices = ntuple(i -> i == d ? range : Colon(), ndims(x))
    return x[indices...]
end
function axis_shift(x, d, shift)
    return circshift(x, ntuple(i -> i == d ? shift : 0, ndims(x)))
end
function upsample_axis(x::AbstractArray{T, N}, d::Int) where {T, N}
    shape = ntuple(i -> i == d ? 1 : size(x, i - (i > d)), N + 1)
    a = reshape(x, shape)

    interleaved = cat(a, zero.(a); dims = d)
    interleaved = x isa Array ? interleaved::Array{T, N + 1} : interleaved

    output_shape = ntuple(
        i -> i == d ? 2 * size(x, i) : size(x, i),
        N,
    )

    return reshape(interleaved, output_shape)
end
function synthesis_axis(x, h, d)
    half = size(x, d) ÷ 2
    low = upsample_axis(axis_slice(x, d, 1:half), d)
    high = upsample_axis(axis_slice(x, d, (half + 1):(2half)), d)
    terms = ntuple(length(h)) do r
        h[r] .* axis_shift(low, d, r - 1) .+
            (isodd(r) ? h[r] : -h[r]) .* axis_shift(high, d, 2 - r)
    end
    return reduce(+, terms)
end
function analysis_axis(x, h, d)
    odd = 1:2:size(x, d)
    lows = ntuple(r -> h[r] .* axis_slice(axis_shift(x, d, 1 - r), d, odd), length(h))
    highs = ntuple(
        r -> (isodd(r) ? h[r] : -h[r]) .* axis_slice(axis_shift(x, d, r - 2), d, odd),
        length(h),
    )
    return cat(reduce(+, lows), reduce(+, highs); dims = d)::typeof(x)
end
function approximation(x)
    return x[ntuple(d -> 1:(size(x, d) ÷ 2), ndims(x))...]
end
function replace_approximation(x::AbstractArray{T, N}, low) where {T, N}
    return _replace_approximation(x, low, Val(1), Val(N))
end
function _replace_approximation(
        x::AbstractArray{T, N},
        low,
        ::Val{D},
        ::Val{N},
    ) where {T, N, D}
    D > N && return low

    n = size(x, D) ÷ 2
    a = axis_slice(x, D, 1:n)
    b = axis_slice(x, D, (n + 1):(2n))

    y = _replace_approximation(a, low, Val(D + 1), Val(N))
    return cat(y, b; dims = D)::typeof(x)
end
function wavelet_synthesis(x, h, levels)
    levels == 0 && return copy(x)
    low = wavelet_synthesis(approximation(x), h, levels - 1)::typeof(x)
    y = replace_approximation(x, low)
    for d in 1:ndims(x)
        y = synthesis_axis(y, h, d)
    end
    return y
end
function wavelet_analysis(x, h, levels)
    levels == 0 && return copy(x)
    y = x
    for d in ndims(x):-1:1
        y = analysis_axis(y, h, d)
    end
    low = wavelet_analysis(approximation(y), h, levels - 1)::typeof(y)
    return replace_approximation(y, low)
end
