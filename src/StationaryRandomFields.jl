module StationaryRandomFields

# Import Modules
using Base
using DocStringExtensions
using FFTW
using Random
using Statistics

# Noise Signal Models
include("./noisesignals/abstractnoisesignal.jl")
include("./noisesignals/noisesignal.jl")
include("./noisesignals/abstractcns.jl")
include("./noisesignals/continuousnoisesignal.jl")

# Power Spectrum Models
include("./powerspectra/abstractpowerspectrum.jl")
include("./powerspectra/modifiers.jl")
include("./powerspectra/singlepowerlaw.jl")
include("./powerspectra/saturatedpowerlaw.jl")
include("./powerspectra/phasepowerlaw.jl")

# Noise Generators
include("./noisegenerators/psnoisegenerator.jl")
include("./noisegenerators/psnoisegen.jl")
end
