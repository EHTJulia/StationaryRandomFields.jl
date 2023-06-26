export ContinuousNoiseSignal

"""
    $(TYPEDEF)

Data type for continuous noise signals. The data type contains an 
AbstractNoiseSignal and defines its corresponding dimensions, RFFT plan, and inverse RFFT plan. 
The type of data (default T=Float64) may also be specified. 

# Fields
$(FIELDS)
"""

struct ContinuousNoiseSignal{T, N <: AbstractNoiseSignal, D <: Tuple, P, PI} <: AbstractContinuousNoiseSignal
    noisesignal::N
    dims::D
    plan::P
    invplan::PI
    function ContinuousNoiseSignal(noisesignal::N, T=Float64) where {N <: AbstractNoiseSignal}
        dims = noisesignal.dims
        temp = zeros(T, dims...)
        plan = plan_rfft(temp)
        invplan = inv(plan)
        return new{T, N, Tuple, typeof(plan), typeof(invplan)}(noisesignal, dims, plan, invplan)
    end
end