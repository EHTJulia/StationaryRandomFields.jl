module PowerSpectrumNoises

# Import Modules
using Base
using FFTW
using Random
using Statistics

# Noise Signal Models
include("./noisesignals/abstractnoisesignal.jl")
include("./noisesignals/noisesignal.jl")

# Power Spectrum Models
include("./powerspectra/abstractpowerspectrum.jl")
include("./powerspectra/singlepowerlaw.jl")
# Noise Generators 
include("noisegenerators/psnoisegenerator.jl")
include("noisegenerators/psnoisegen1d.jl")
include("noisegenerators/psnoisegen2d.jl")
end
