export NoiseSignal

"""
    $(TYPEDEF)

This is a data type for noise signals of any dimension.
"""
struct NoiseSignal <: AbstractNoiseSignal
    dims::Tuple
end