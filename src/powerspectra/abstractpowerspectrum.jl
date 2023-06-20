export AbstractPowerSpectrumModel

"""
    AbstractPowerSpectrumModel{N}

The abstract model type for a PowerSpectrum at `N` dimension. To instantiate your own model type you should
subtybe from this model. Additionally you need to implement the following
methods to satify the interface:

**Mandatory Methods**
- [`isprimitive`](@ref): defines whether a model is standalone or is defined in terms of other models.
   is the model is primitive then this should return `IsPrimitive()` otherwise it returns
   `NotPrimitive()`
- [`fourieranalytic`](@ref): defines whether the model fourier coefficients  can be computed analytically. If yes
   then this should return `IsAnalytic()` and the user *must* to define `visibility_point`.
   If not analytic then `visanalytic` should return `NotAnalytic()`.
- [`fourier_point`](@ref): Defines how to compute model fourier coefficients pointwise. Note this is
    must be defined if `visanalytic(::Type{YourModel})==IsAnalytic()`.
"""
abstract type AbstractPowerSpectrumModel{N} end

# return the number of dimensions
Base.ndims(::AbstractPowerSpectrumModel{N}) where {N} = N

"""
    DensityAnalytic

Internal type for specifying the nature of the model functions.
Whether they can be easily evaluated pointwise analytic. This
is an internal type that may change.
"""
abstract type DensityAnalytic end

"""
    $(TYPEDEF)

Defines a trait that a states that a model is analytic.
This is usually used with an abstract model where we use
it to specify whether a model has a analytic fourier transform
and/or image.
"""
struct IsAnalytic <: DensityAnalytic end

"""
    $(TYPEDEF)

Defines a trait that a states that a model is analytic.
This is usually used with an abstract model where we use
it to specify whether a model has does not have a easy analytic
fourier transform and/or intensity function.
"""
struct NotAnalytic <: DensityAnalytic end

"""
    fourieranalytic(::Type{<:AbstractPowerSpectrumModel})

Determines whether the model is pointwise analytic in Fourier domain, i.e. we can evaluate
its fourier transform at an arbritrary point.

If `IsAnalytic()` then it will try to call `visibility_point` to calculate the complex visibilities.
Otherwise it fallback to using the FFT that works for all models that can compute an image.

"""
@inline fourieranalytic(::AbstractPowerSpectrumModel) = NotAnalytic()

#=
    This is internal function definitions for how to
    compose whether a model is analytic. We need this
    for composite models.
=#
Base.@constprop :aggressive Base.:*(::IsAnalytic, ::IsAnalytic) = IsAnalytic()
Base.@constprop :aggressive Base.:*(::IsAnalytic, ::NotAnalytic) = NotAnalytic()
Base.@constprop :aggressive Base.:*(::NotAnalytic, ::IsAnalytic) = NotAnalytic()
Base.@constprop :aggressive Base.:*(::NotAnalytic, ::NotAnalytic) = NotAnalytic()

"""
    fourier_point(model::AbstractPowerSpectrumModel, p)

Function that computes the pointwise fourier component.
This must be implemented in the model interface
if `fourieranalytic(::Type{MyModel}) == IsAnalytic()`
"""
function fourier_point end
