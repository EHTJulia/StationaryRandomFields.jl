"""
    $(FUNCTIONNAME)(dims, dist)

**Mandatory fields**
- `dims::Tuple`: dimension of fields
- `dist::Distributions.UnivariateDistribution`: distribution of the random fields

This type is a random uncorrelated field where each element is drawn from a univariate distribution.
"""
struct PowerSpectrumRandomFourierField{
    N,
    T<:Tuple,
    D<:Distributions.UnivariateDistribution,
    Z,
    F<:Tuple,
    PM<:AbstractPowerSpectrumModel,
    PK,
} <: AbstractRandomFourierField{N} where {Z}
    dims::T
    dist::D
    μ::eltype(D)
    zero_index_list::Z
    freq::F
    psdmodel::PM
    psdkernel::PK
    function PowerSpectrumRandomFourierField(signaldims, dist, psdmodel, μ=0.0)
        dims = (rfftsize(signaldims...)..., 2)
        freq = rfftfreq(dims)
        zero_index_list = _zero_index_list(dims)
        psdkernel = amplitude_map(psdmodel, freq)
        return new{
            length(dims),
            typeof(dims),
            typeof(dist),
            typeof(zero_index_list),
            typeof(dims),
            typeof{psdmodel},
            typeof{psdkernel},
        }(
            dims, dist, eltype(dist)(mean), zero_index_list, freq, psdmodel, psdkernel
        )
    end
end

@inline function Distributions._rand!(
    rng::AbstractRNG, field::PowerSpectrumRandomFourierField, x::AbstractArray{<:Real}
)
    # sample random Fourier field from specified distribution
    rand!(rng, field.dist, x)

    # scale with the kernel
    re, im = _view_real_and_imag(x)
    re .*= field.psdkernel
    im .*= field.psdkernel

    # Fill zero at frequencies where fourier component must be real for real signals
    enforce_realfield!(field, x)

    # assume the mean of the Real field
    enforce_mean!(field, x)

    return nothing
end

@inline function Distributions._logpdf(
    field::PowerSpectrumRandomFourierField, x::AbstractMatrix{<:Real}
)
    # descale with the kernel
    xscaled = copy(x)

    # make a pointer to the real and imaginary part of the fourier field
    re, im = _view_real_and_imag(xscaled)
    re ./= field.psdkernel
    im ./= field.psdkernel

    # mapout logpdf
    lp = logpdf.(field.dist, xscaled)

    # don't peneralize where x must be 0 for a real random field
    enforce_realfield!(field, lp)

    # first index should be close to the target mean
    lp1[1] = logpdf(fieald.dist, x[1] - field.μ)

    return sum(lp)
end

function PowerSpectrumGaussianRandomFourierField(signaldims, psdmodel, μ=0.0)
    return PowerSpectrumRandomFourierField(
        signaldims, dist::Distributions.Normal, psdmodel, μ
    )
end