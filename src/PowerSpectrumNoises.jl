module PowerSpectrumNoises

# Import Modules
using Base
using DocStringExtensions
using FFTW
using Random
using Statistics

# Noise Signal Models
include("./noisesignals/abstractnoisesignal.jl")
include("./noisesignals/noisesignal1d.jl")
include("./noisesignals/noisesignal2d.jl")

# Power Spectrum Models
include("./powerspectra/abstractpowerspectrum.jl")
include("./powerspectra/modifiers.jl")

# Noise Generators
#include("noisegenerators/psnoisegenerator.jl")
#include("noisegenerators/psnoisegen1d.jl")
#include("noisegenerators/psnoisegen2d.jl")
end
