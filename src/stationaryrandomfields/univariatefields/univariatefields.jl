export UnivariateRandomField
export UnivariateGaussianRandomField
export UnivariateLaplaceRandomField

"""
    $(FUNCTIONNAME)(dims, dist)

**Mandatory fields**
- `dims::Tuple`: dimension of fields
- `dist::Distributions.UnivariateDistribution`: distribution of the random fields

This type is a random uncorrelated field where each element is drawn from a univariate distribution.
"""
struct UnivariateRandomField{N,T<:Tuple,D<:Distributions.UnivariateDistribution} <:
       AbstractUnivariateRandomField{N}
    dims::T
    dist::D
    function UnivariateRandomField(dims, dist)
        return new{length(dims),typeof(dims),typeof(dist)}(dims, dist)
    end
end

@inline Base.eltype(field::UnivariateRandomField) = eltype(field.dist)

"""
    $(TYPEDEF)

This type is a random field where each element is drawn from a univariate Gaussian (i.e. Normal) distribution.
"""
@inline UnivariateGaussianRandomField(
    dims, dist::Distributions.Normal=DefaultGaussianDistribution
) = UnivariateRandomField(dims, dist)

"""
    $(TYPEDEF)

This type is a random field where each element is drawn from a univariate Gaussian (i.e. Normal) distribution.
"""
@inline UnivariateLaplaceRandomField(
    dims, dist::Distributions.Laplace=DefaultLaplaceDistribution
) = UnivariateRandomField(dims, dist)