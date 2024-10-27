export AbstractRandomFourierField

"""
    $(TYPEDEF)

This is an abstract field type for the `N`-dimensional stationary uncorrelated random fields where each element is drawn from a univariate distribution.
This particular abstract type assumes that the field is given by the Fourier transform of a real field on the Signal domain.
Given that such a field is described by Hermite-conjugate complex numbers, `dims` should be a tuple of integers with the last element being 2 (real and imag).

The field type contains the information for the dimension of the field and the distribution of the random fields.

**Mandatory fields**
- `dims::Tuple`: dimension of fields. The last dimension should be 2, corresponding to the real and imaginary parts of the complex Fourier field.
- `dist::Distributions.UnivariateDistribution`: distribution of the random fields in Fourier Domain
- `μ::eltype(dist)`: mean of the field in the Signal domain
- `zero_index_list`: list of indices where the Fourier field must be real numbers`
- `freq`: list of frequencies in each of the dimensions`

**Mandatory methods**
See the documentation of `AbstractUnivariateRandomField` for the methods required as a subtype of AbstractStationaryRandomField.
"""
abstract type AbstractRandomFourierField{N} <: AbstractStationaryRandomField{N} end


@inline function enforce_mean!(field::AbstractRandomFourierField, x::AbstractArray{<:Real})
    x[1] = field.μ
    return nothing
end


@inline function enforce_realfield!(field::AbstractRandomFourierField, x::AbstractArray{<:Real})
    @assert size(x) == field.dims

    # get a view on the imaginary part
    ϵim = _view_imag(x)

    # zero-frequency point must be real
    zero_const = zero(eltype(field))
    for idx in fields.zero_index_list
        ϵim[idx...] = zero_const
    end

    return nothing
end

# 
# Helper functions
#
"""
    $(FUNCTIONNAME)(dims::Tuple)


Find the frequency index where fourier field signal must be real number to make the real field
"""
@inline function _zero_index_list(dims::Tuple)
    index_list = []
    
    push!(index_list, (1,))

    evenness = iseven.(dims)
    midindex = Tuple(dim ÷ 2 + 1 for dim in dims)
    if evenness[1]
        # create iterator
        iter = []
        for idim in 1:length(dims)-1
            if evenness[idim]
                push!(iter, tuple(1, midindex[idim]))
            else
                push!(iter, tuple(1,))
            end
        end

        # fill zero
        for idx in Iterators.product(iter...)
            push!(index_list, idx)
        end
    end
    return index_list
end


@inline function _view_real_and_imag(x::AbstractArray{<:Real})
    viewidx = []
    for idim in 1:ndims(x)-1
        push!(viewidx, :)
    end
    ϵre = view(x, viewidx..., 1)
    ϵim = view(x, viewidx..., 2)
    return ϵre, ϵim
end


@inline function _view_real(x::AbstractArray{<:Real})
    viewidx = []
    for idim in 1:ndims(x)-1
        push!(viewidx, :)
    end
    ϵre = view(x, viewidx..., 1)
    return ϵre
end


@inline function _view_imag(x::AbstractArray{<:Real})
    viewidx = []
    for idim in 1:ndims(x)-1
        push!(viewidx, :)
    end
    ϵim = view(x, viewidx..., 2)
    return ϵim
end