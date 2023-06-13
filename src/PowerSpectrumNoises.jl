module PowerSpectrumNoises

# Import Modules
using Base
using FFTW
using Random
using Statistics

# Noise Signal Models
include("./noisesignals/noisesignal1d.jl")
include("./noisesignals/noisesignal2d.jl")
include("./noisesignals/abstractnoisesignal.jl")

# Power Spectrum Models
include("./powerspectra/abstractpowerspectrum.jl")
include("./powerspectra/singlepowerlaw1d.jl")

# Noise Generators 
include("noisegenerators/psnoisegenerator.jl")
include("noisegenerators/psnoisegen1d.jl")
end
