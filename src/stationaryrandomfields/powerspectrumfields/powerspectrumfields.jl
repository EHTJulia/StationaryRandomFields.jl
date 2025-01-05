"""
    $(FUNCTIONNAME)(dims, dist)

**Mandatory fields**
- `dims::Tuple`: dimension of fields
- `dist::Distributions.UnivariateDistribution`: distribution of the random fields

This type is a random uncorrelated field where each element is drawn from a univariate distribution.
"""
struct PowerSpectrumRandomField{N,T<:Tuple,D<:AbstractRandomFourierField,P,PI} <:
       AbstractPowerSpectrumRandomField{N}
    dims::T
    dist::D
    fftplan::P
    ifftplan::PI
    function PowerSpectrumRandomField(dims, dist, psdmodel, μ=0.0)
        # initialize random Fourier field
        fourierfield = PowerSpectrumRandomFourierField(dims, dist, psdmodel, μ)

        # initialize fft and ifft plans
        T = eltype(dist)
        fftplan = plan_rfft(Array{T}(undef, dims...))
        ifftplan = inv(fftplan)

        return new{
            length(dims),typeof(dims),typeof(fourierfield),typeof(fftplan),typeof(ifftplan)
        }(
            dims, fourierfield, fftplan, ifftplan
        )
    end
end

function PowerSpectrumGaussianRandomField(dims, psdmodel, μ=0.0)
    return PowerSpectrumRandomField(dims, dist::Distributions.Normal, psdmodel, μ)
end