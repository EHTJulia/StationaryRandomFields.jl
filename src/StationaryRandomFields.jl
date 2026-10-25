module StationaryRandomFields

using Distributions
using FastHartleyTransform
using LinearAlgebra
using Random
using Statistics
using Wavelets

import Distributions: _rand!, _logpdf, insupport, logpdf

export AbstractStationaryRandomField,
    UnivariateRandomField,
    UnivariateGaussianRandomField,
    UnivariateLaplaceRandomField,
    TransformedRandomField,
    AbstractFieldTransform,
    IdentityTransform,
    HartleyTransform,
    SpectralScaling,
    DynamicSpectralScaling,
    WaveletTransform,
    ScaleTransform,
    TransformChain,
    FieldGrid,
    prepare,
    apply,
    apply!,
    input_size,
    output_size,
    AbstractPowerSpectrumModel,
    SinglePowerLaw,
    MaternPS,
    MarkovPS,
    SaturatedPowerLaw,
    ConstantPowerSpectrum,
    power_point,
    amplitude_point,
    power_map,
    amplitude_map,
    power_map!,
    amplitude_map!,
    stretched_ps,
    rotated_ps,
    renormed_ps,
    modify_ps,
    PowerSpectrumStretch,
    PowerSpectrumRotate,
    PowerSpectrumRenormalize,
    PeriodicCovariance,
    covariance_value,
    spectral_kernel

include("fields.jl")
include("spectra.jl")
include("transforms.jl")
include("wavelets.jl")
include("transformed.jl")
include("covariance.jl")
include("display.jl")
include("descriptors.jl")

end
