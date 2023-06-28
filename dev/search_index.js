var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = PowerSpectrumNoises","category":"page"},{"location":"#PowerSpectrumNoises","page":"Home","title":"PowerSpectrumNoises","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for PowerSpectrumNoises.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [PowerSpectrumNoises]","category":"page"},{"location":"#PowerSpectrumNoises.AbstractNoiseSignal","page":"Home","title":"PowerSpectrumNoises.AbstractNoiseSignal","text":"AbstractNoiseSignal\n\nThis is an abstract data type for signal data on which powerspectrum noises will be created. The data type contains information to define the dimensions and sizes of data, and provides various methods.\n\nMandatory fields\n\ndims::Tuple: dimension of data\n\nMandatory methods\n\nsize: returns the size of the signal\nsizeof: returns product of signal dimensions\nndims: returns the number of the dimension\nrfftsize: returns the size of the signal in the Fourier plane (for rfft)\nrfftfreq: returns the frequency grid along each dimension in the Fourier plane (for rfft)\ngenerate_gaussian_noise: returns Gaussian noises in Fourier domain with the size of rfftsize(data)\n\n\n\n\n\n","category":"type"},{"location":"#PowerSpectrumNoises.AbstractPowerSpectrumModel","page":"Home","title":"PowerSpectrumNoises.AbstractPowerSpectrumModel","text":"AbstractPowerSpectrumModel{N}\n\nThe abstract model type for a PowerSpectrum at N dimension. To instantiate your own model type you should subtype from this model. Additionally you need to implement the following methods to satify the interface:\n\nMandatory Methods\n\nfourieranalytic: defines whether the model fourier coefficients  can be computed analytically. If yes  then this should return IsAnalytic() and the user must to define visibility_point.  If not analytic then fourieranalytic should return NotAnalytic().\npower_point: Defines how to compute model power of fourier coefficients pointwise. Note this is   must be defined if fourieranalytic(::Type{YourModel})==IsAnalytic().\n\nOptional Methods\n\namplitude_point: Defines how to compute model amplitudes of fourier coefficients pointwise. Note this is   must be defined if fourieranalytic(::Type{YourModel})==IsAnalytic().   It should be defined as √(power_point(model, ν...)).\n\n\n\n\n\n","category":"type"},{"location":"#PowerSpectrumNoises.DensityAnalytic","page":"Home","title":"PowerSpectrumNoises.DensityAnalytic","text":"DensityAnalytic\n\nInternal type for specifying the nature of the model functions. Whether they can be easily evaluated pointwise analytic. This is an internal type that may change.\n\n\n\n\n\n","category":"type"},{"location":"#PowerSpectrumNoises.IsAnalytic","page":"Home","title":"PowerSpectrumNoises.IsAnalytic","text":"struct IsAnalytic <: PowerSpectrumNoises.DensityAnalytic\n\nDefines a trait that a states that a model is analytic. This is usually used with an abstract model where we use it to specify whether a model has a analytic fourier transform and/or image.\n\n\n\n\n\n","category":"type"},{"location":"#PowerSpectrumNoises.ModelModifier","page":"Home","title":"PowerSpectrumNoises.ModelModifier","text":"abstract type ModelModifier{T}\n\nGeneral type for a model modifier. These transform any model using simple Fourier transform properties. To modify a model you can use the ModifiedPowerSpectrumModel{N} constructor or the modify function.\n\njulia> visanalytic(stretched(Disk(), 2.0, 2.0)) == visanalytic(Disk())\ntrue\n\nTo implement a model transform you need to specify the following methods:\n\ntransform_ν\ntransform_image\nscale_fourier\nscale_image\nradialextent\n\nSee the docstrings of those methods for guidance on implementation details.\n\nAdditionally these methods assume the modifiers are of the form\n\nI(x,y) -> fᵢ(x,y)I(gᵢ(x,y)) V(u,v) -> fᵥ(u,v)V(gᵥ(u,v))\n\nwhere g are the transformimage/uv functions and f are the scaleimage/uv function.\n\n\n\n\n\n","category":"type"},{"location":"#PowerSpectrumNoises.ModifiedPowerSpectrumModel","page":"Home","title":"PowerSpectrumNoises.ModifiedPowerSpectrumModel","text":"struct ModifiedPowerSpectrumModel{N, M<:AbstractPowerSpectrumModel{N}, T<:Tuple} <: AbstractPowerSpectrumModel{N}\n\nContainer type for models that have been transformed in some way. For a list of potential modifiers or transforms see subtypes(ModelModifiers).\n\nFields\n\nmodel: base model\ntransform: model transforms\n\n\n\n\n\n","category":"type"},{"location":"#PowerSpectrumNoises.NoiseSignal","page":"Home","title":"PowerSpectrumNoises.NoiseSignal","text":"struct NoiseSignal{T<:Tuple} <: AbstractNoiseSignal\n\nThis is a data type for noise signals of any dimension.\n\n\n\n\n\n","category":"type"},{"location":"#PowerSpectrumNoises.NotAnalytic","page":"Home","title":"PowerSpectrumNoises.NotAnalytic","text":"struct NotAnalytic <: PowerSpectrumNoises.DensityAnalytic\n\nDefines a trait that a states that a model is analytic. This is usually used with an abstract model where we use it to specify whether a model has does not have a easy analytic fourier transform and/or intensity function.\n\n\n\n\n\n","category":"type"},{"location":"#PowerSpectrumNoises.Renormalize","page":"Home","title":"PowerSpectrumNoises.Renormalize","text":"struct Renormalize{T} <: PowerSpectrumNoises.ModelModifier{T}\n\nRenormalizes the flux of the model to the new value scale*flux(model). We have also overloaded the Base.:* operator as syntactic sugar although I may get rid of this.\n\nExample\n\njulia> modify(Gaussian(), Renormalize(2.0)) == 2.0*Gaussian()\ntrue\n\n\n\n\n\n","category":"type"},{"location":"#PowerSpectrumNoises.Rotate","page":"Home","title":"PowerSpectrumNoises.Rotate","text":"Rotate(ξ)\n\nType for the rotated model. This is more fine grained control of rotated model.\n\nExample\n\njulia> modify(Gaussian(), Rotate(2.0)) == rotated(Gaussian(), 2.0)\ntrue\n\n\n\n\n\n","category":"type"},{"location":"#PowerSpectrumNoises.SaturatedPowerLaw","page":"Home","title":"PowerSpectrumNoises.SaturatedPowerLaw","text":"struct SaturatedPowerLaw{N, T} <: AbstractPowerSpectrumModel{N}\n\nType for power law with input dimension {N} and profile\n\n    P(\nu) =  \nu  ^eta ndcases\n\nwhen inscale < |ν| < outscale. β refers to the input index. \n\nWhen |ν| > outscale: P(ν) = 0, and when |ν| < inscale:\n\n    P(\nu) = (inscale) ^eta ndcases\n\nCan be renormalized, stretched, and rotated via ModelModifier.\n\nFields\n\nindex\ninscale\noutscale\n\n\n\n\n\n","category":"type"},{"location":"#PowerSpectrumNoises.SinglePowerLaw","page":"Home","title":"PowerSpectrumNoises.SinglePowerLaw","text":"struct SinglePowerLaw{N, T} <: AbstractPowerSpectrumModel{N}\n\nType for single power law with input dimension {N} and profile\n\n    P(\nu) =  \nu  ^eta ndcases\n\nwhere β is the input index and ν is frequency. Can be renormalized, stretched, and rotated via ModelModifier.\n\nFields\n\nindex\n\n\n\n\n\n","category":"type"},{"location":"#PowerSpectrumNoises.Stretch","page":"Home","title":"PowerSpectrumNoises.Stretch","text":"Stretch(α, β)\nStretch(r)\n\nStretched the model in the x and y directions, i.e. the new intensity is     Iₛ(x,y) = 1/(αβ) I(x/α, y/β), where were renormalize the intensity to preserve the models flux.\n\nExample\n\njulia> modify(Gaussian(), Stretch(2.0)) == stretched(Gaussian(), 2.0, 1.0)\ntrue\n\nIf only a single argument is given it assumes the same stretch is applied in both direction.\n\njulia> Stretch(2.0) == Stretch(2.0, 2.0)\ntrue\n\n\n\n\n\n","category":"type"},{"location":"#PowerSpectrumNoises.amplitude_map-Tuple{AbstractPowerSpectrumModel, Any}","page":"Home","title":"PowerSpectrumNoises.amplitude_map","text":"amplitude_map(model::AbstractPowerSpectrumModel, data)\n\nFunction that maps the amplitude function of a signal data frequency grid.  The intended input \"data\" is the frequency grid output by FFTW.rfftfreq.  Alternatively a SignalNoise or ContinuousSignalNoise object can be input, and the corresponding frequency grid will be computed and mapped.\n\nIn default, it should be defined as .√(power_point(model, data)).\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.amplitude_point-Tuple{AbstractPowerSpectrumModel, Vararg{Any}}","page":"Home","title":"PowerSpectrumNoises.amplitude_point","text":"amplitude_point(model::AbstractPowerSpectrumModel, ν...)\n\nFunction that computes the pointwise amplitude of fourier component at a set of frequencies ν. This must be implemented in the model interface if fourieranalytic(::Type{MyModel}) == IsAnalytic().\n\nIn default, it should be defined as √(power_point(model, ν...)).\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.basemodel-Tuple{PowerSpectrumNoises.ModifiedPowerSpectrumModel}","page":"Home","title":"PowerSpectrumNoises.basemodel","text":"basemodel(model::ModifiedPowerSpectrumModel)\n\nReturns the ModifiedPowerSpectrumModel with the last transformation stripped.\n\nExample\n\njulia> basemodel(stretched(Disk(), 1.0, 2.0)) == Disk()\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.fourieranalytic-Tuple{AbstractPowerSpectrumModel}","page":"Home","title":"PowerSpectrumNoises.fourieranalytic","text":"fourieranalytic(::Type{<:AbstractPowerSpectrumModel})\n\nDetermines whether the model is pointwise analytic in Fourier domain, i.e. we can evaluate its fourier transform at an arbritrary point.\n\nIf IsAnalytic() then it will try to call visibility_point to calculate the complex visibilities. Otherwise it fallback to using the FFT that works for all models that can compute an image.\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.generate_gaussian_noise-Tuple{AbstractNoiseSignal}","page":"Home","title":"PowerSpectrumNoises.generate_gaussian_noise","text":"generate_gaussian_noise(data::AbstractNoiseSignal; rng=Random.default_rng())\n\nGenerate complex gausisian noises on Fourier space over the frequency space expected for Real FFT. This return a complex array, which is consistent with Fourier-domain representation of a real uncorrelated Gaussian noise with the zero mean and variance of sizeof(data)/2 in signal-domain. The size of the output array is given by rfftsize(data).\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.get_fourier_noise-Tuple{PowerSpectrumNoiseGenerator, AbstractArray}","page":"Home","title":"PowerSpectrumNoises.get_fourier_noise","text":"get_fourier_noise(psgen::PowerSpectrumNoiseGenerator, signoise::AbstractArray)\n\nRetrieve power-spectrum-scaled fourier noise from given signal noise by transforming to fourier space.\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.get_power_spectrum-Tuple{PowerSpectrumNoiseGenerator, AbstractArray}","page":"Home","title":"PowerSpectrumNoises.get_power_spectrum","text":"get_power_spectrum(psgen::PowerSpectrumNoiseGenerator, signoise::AbstractArray)\n\nRetrieve power spectrum from given signal noise.\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.modify-Tuple{AbstractPowerSpectrumModel, Vararg{Any}}","page":"Home","title":"PowerSpectrumNoises.modify","text":"modify(m::AbstractPowerSpectrumModel, transforms...)\n\nModify a given model using the set of transforms. This is the most general function that allows you to apply a sequence of model transformation for example\n\nmodify(Gaussian(), Stretch(2.0, 1.0), Rotate(π/4), Shift(1.0, 2.0), Renorm(2.0))\n\nwill create a asymmetric Gaussian with position angle π/4 shifted to the position (1.0, 2.0) with a flux of 2 Jy. This is similar to Flux's chain.\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.posangle-Tuple{Rotate}","page":"Home","title":"PowerSpectrumNoises.posangle","text":"posangle(model)\n\n\nReturns the rotation angle of the rotated model\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.power_map-Tuple{AbstractPowerSpectrumModel, Tuple}","page":"Home","title":"PowerSpectrumNoises.power_map","text":"power_map(model, gridofν)\n\nFunction that maps the power law function of a signal data frequency grid.  The intended input \"data\" is the frequency grid output by FFTW.rfftfreq.  Alternatively a SignalNoise or ContinuousSignalNoise object can be input,  and the corresponding frequency grid will be computed and mapped.\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.power_point","page":"Home","title":"PowerSpectrumNoises.power_point","text":"power_point(model::AbstractPowerSpectrumModel, ν...)\n\nFunction that computes the pointwise amplitude of fourier component at a set of frequencies ν. This must be implemented in the model interface if fourieranalytic(::Type{MyModel}) == IsAnalytic()\n\n\n\n\n\n","category":"function"},{"location":"#PowerSpectrumNoises.renormed-Union{Tuple{M}, Tuple{M, Any}} where M<:AbstractPowerSpectrumModel","page":"Home","title":"PowerSpectrumNoises.renormed","text":"renormed(model, f)\n\n\nRenormalizes the model m to have total flux f*flux(m). This can also be done directly by calling Base.:* i.e.,\n\njulia> renormed(m, f) == f*M\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.rfftsize-Tuple{AbstractNoiseSignal}","page":"Home","title":"PowerSpectrumNoises.rfftsize","text":"rfftsize(data::AbstractNoiseSignal)::Tuple\n\nReturns the size of the Fourier-domain signal\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.rotated-Tuple{Any, Any}","page":"Home","title":"PowerSpectrumNoises.rotated","text":"rotated(model, ξ)\n\n\nRotates the model by an amount ξ in radians in the clockwise direction.\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.scale_fourier","page":"Home","title":"PowerSpectrumNoises.scale_fourier","text":"scale_fourier(model::AbstractModifier, ν...)\n\nReturns a number on how to scale the fourier coefficients at the frequency coordinate ν for an modified model\n\n\n\n\n\n","category":"function"},{"location":"#PowerSpectrumNoises.stretched-Tuple{AbstractPowerSpectrumModel, Vararg{Any}}","page":"Home","title":"PowerSpectrumNoises.stretched","text":"stretched(model, α)\n\n\nStretches the model m according to the formula     Vₛ(u) = V(u1/α1, u2/α2, ...), where were renormalize the intensity to preserve the models flux.\n\n\n\n\n\n","category":"method"},{"location":"#PowerSpectrumNoises.transform_ν","page":"Home","title":"PowerSpectrumNoises.transform_ν","text":"transform_ν(model::AbstractModifier, ν...)\n\nReturns a transformed frequency coordinate ν according to the model modifier\n\n\n\n\n\n","category":"function"},{"location":"#PowerSpectrumNoises.unmodified-Tuple{PowerSpectrumNoises.ModifiedPowerSpectrumModel}","page":"Home","title":"PowerSpectrumNoises.unmodified","text":"unmodified(model::ModifiedPowerSpectrumModel)\n\nReturns the un-modified model\n\nExample\n\njulia> m = stretched(rotated(Gaussian(), π/4), 2.0, 1.0)\njulia> umodified(m) == Gaussian()\ntrue\n\n\n\n\n\n","category":"method"}]
}
