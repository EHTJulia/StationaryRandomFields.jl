module StationaryRandomFields

# Import Modules
using Base
using Distributions
import Distributions
using DocStringExtensions
using FFTW
using Random
using Statistics

# Set the default noise sampler in Fourier space
const DefaultFourierDistribution = Distributions.Normal(0, 1/√2)
const DefaultGaussianDistribution = Distributions.Normal()
const DefaultLaplaceDistribution = Distributions.Laplace()

# AbstractStationaryRandomField
include("stationaryrandomfields/abstract.jl")

# Power Spectrum Models
include("./powerspectra/abstractpowerspectrum.jl")
include("./powerspectra/modifiers.jl")
include("./powerspectra/singlepowerlaw.jl")
include("./powerspectra/saturatedpowerlaw.jl")

# Univariate Random Fields
include("stationaryrandomfields/univariatefields/abstract.jl")
include("stationaryrandomfields/univariatefields/univariatefields.jl")

# Univariate Random Fourier Fields
include("stationaryrandomfields/fourierfields/abstract.jl")
include("stationaryrandomfields/fourierfields/rfftutils.jl")
include("stationaryrandomfields/fourierfields/univariatefields.jl")
include("stationaryrandomfields/fourierfields/powerspectrumfields.jl")

# Powerspectrum Random Fourier Fields
include("stationaryrandomfields/powerspectrumfields/abstract.jl")
include("stationaryrandomfields/powerspectrumfields/powerspectrumfields.jl")

# Noise Signal Models
#include("./noisesignals/abstractnoisesignal.jl")
#include("./noisesignals/noisesignal.jl")
#include("./noisesignals/abstractcns.jl")
#include("./noisesignals/continuousnoisesignal.jl")

# Noise Generators
#include("./noisegenerators/psnoisegenerator.jl")
#include("./noisegenerators/psnoisegen.jl")
end
