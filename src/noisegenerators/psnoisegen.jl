export PSNoiseGenerator

"""
    $(TYPEDEF)

Data type which takes a power spectrum model (psmodel) and continuous noise signal (noisesignal) 
and generates correlated signal noise.

# Fields
$(FIELDS)
"""

struct PSNoiseGenerator{P<:AbstractPowerSpectrumModel, N<:ContinuousNoiseSignal} <: PowerSpectrumNoiseGenerator
    psmodel::P
    noisesignal::N
end

