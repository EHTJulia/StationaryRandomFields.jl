export AbstractContinuousNoiseSignal

"""
AbstractContinuousNoiseSignal

This is an abstract data type for N dimensional continuous signal data on which signal noise will be generated. 
The data type contains an AbstractNoiseSignal and defines its corresponding dimensions, RFFT plan, and inverse RFFT plan. 
The type of data (default Float64) may also be specified.

"""

abstract type AbstractContinuousNoiseSignal end

