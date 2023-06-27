export NoiseSignal

"""
    $(TYPEDEF)

This is a data type for noise signals of any dimension.
"""
struct NoiseSignal{T <: Tuple} <: AbstractNoiseSignal
    dims::T
end