# Introduction

`StationaryRandomFields.jl` simulates realistic correlated noise for signal data of any given dimensions. 
The package follows the power-law noise procedure introduced by [Timmer et al. (1995)](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..707T/abstract); random Gaussian noise is drawn in the Fourier frequency domain and scaled by the square root of a power law spectrum. 
An inverse transform back to the signal domain gives the stochastic power-law noise. 

The module currently
- Provides abstract data types and methods to define noise signals, construct and modify power-law scaling functions, and generate signal noise
- Implements multiple power spectrum types:
    - Basic spectra of form $P = \nu^{\alpha}$ ([`SinglePowerLaw`](@ref))
    - Cut-off spectra with inner and outer scales ([`SaturatedPowerLaw`](@ref))
- Provides methods to reverse the process and retrieve underlying power spectra from input signal noise

## Installation
The package can be installed by running:
```julia
using Pkg
Pkg.add("StationaryRandomFields")
```
