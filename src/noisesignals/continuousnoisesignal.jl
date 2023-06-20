export ContinuousNoiseSignal

"""
    ContinuousNoiseSignal <: AbstractContinuousNoiseSignal

This is a data type for N dimensional continuous noise signals.
"""

struct ContinuousNoiseSignal{T, N <: AbstractNoiseSignal, D <: Tuple, P, PI} <: AbstractContinuousNoiseSignal
    noisesignal::N
    dims::D
    plan::P
    invplan::PI
    function ContinuousNoiseSignal(noisesignal, T = Float64) 
        dims = noisesignal.dims
        temp = zeros(T, dims...)
        plan = plan_rfft(temp)
        invplan = inv(plan)
        return new{T, AbstractNoiseSignal, Tuple, typeof(plan), typeof(invplan)}(noisesignal, dims, plan, invplan)
    end
end