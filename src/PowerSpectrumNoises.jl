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

# Power Spectrum Models
include("./powerspectra/abstractpowerspectrum.jl")
include("./powerspectra/singlepowerlaw1d.jl")

end
