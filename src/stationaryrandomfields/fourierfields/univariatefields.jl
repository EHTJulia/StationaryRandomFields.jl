export UnivariateRandomFourierField
export UnivariateGaussianRandomFourierField

"""
    $(FUNCTIONNAME)(signaldims, dist, μ=0.)

- `signaldims::Tuple`: dimension of the random fields on the signal domain
- `dist::Distributions.UnivariateDistribution`: distribution of the random fields on Fourier domain
- `μ::eltype(dist)`: mean of the field in the Signal domain

This type is a random uncorrelated field where each element is drawn from a univariate distribution.
"""
struct UnivariateRandomFourierField{N, T<:Tuple, D<:Distributions.UnivariateDistribution, Z, F<:Tuple} <: AbstractRandomFourierField{N} where Z
    dims::T
    dist::D
    μ::eltype(D)
    zero_index_list::Z
    freq::F
    function UnivariateRandomFourierField(signaldims, dist, μ=0.)
        dims=(rfftsize(signaldims...)..., 2)
        freq=rfftfreq(dims)
        zero_index_list = _zero_index_list(dims)
        return new{length(dims), typeof(dims), typeof(dist), typeof(zero_index_list), typeof(dims)}(dims, dist, eltype(dist)(mean), zero_index_list, freq)
    end
end


@inline function Distributions._rand!(rng::AbstractRNG, field::UnivariateRandomFourierField, x::AbstractArray{<:Real})
    # sample random Fourier field from specified distribution
    rand!(rng, field.dist, x)

    # assume the mean of the Real field
    enforce_mean!(field, x)

    # Fill zero at frequencies where fourier component must be real for real signals
    enforce_realfield!(field, x)

    return nothing
end


@inline function Distributions._logpdf(field::UnivariateRandomFourierField, x::AbstractMatrix{<:Real})
    # mapout logpdf
    lp = logpdf.(field.dist, x)

    # don't peneralize where x must be 0 for a real random field
    enforce_realfield!(field, lp)

    # first index should be close to the target mean
    lp1[1] = logpdf(fieald.dist, x[1]-field.μ)

    return sum(lp)
end

UnivariateGaussianRandomFourierField(signaldims, μ=0.) = UnivariateRandomFourierField(signaldims, dist::Distributions.Normal, μ)