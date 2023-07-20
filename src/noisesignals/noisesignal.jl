export NoiseSignal

"""
    $(TYPEDEF)

This is a data type for noise signals of any dimension.
    
# Fields
$(FIELDS)
"""
struct NoiseSignal{T <: Tuple} <: AbstractNoiseSignal
    dims::T
end