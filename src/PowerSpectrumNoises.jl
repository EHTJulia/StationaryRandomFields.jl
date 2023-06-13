module PowerSpectrumNoises

# Import Modules
using Base
using FFTW
using Random
using Statistics

# Noise Signal Models
include("./noisesignals/abstractnoisesignal.jl")
include("./noisesignals/noisesignal1d.jl")
include("./noisesignals/noisesignal2d.jl")

# Pusing ower Spectrum Models
include("./powerspectra/abstractpowerspectrum.jl")
include("./powerspectra/singlepowerlaw1d.jl")

# Noise Generators 
include("noisegenerators/psnoisegenerator.jl")
include("noisegenerators/psnoisegen1d.jl")
end
