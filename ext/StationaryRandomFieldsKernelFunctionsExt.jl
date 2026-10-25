module StationaryRandomFieldsKernelFunctionsExt

using KernelFunctions: KernelFunctions
using StationaryRandomFields

struct PeriodicSpectralKernel{C <: PeriodicCovariance} <: KernelFunctions.Kernel
    covariance::C
end
(k::PeriodicSpectralKernel)(x, y) = covariance_value(k.covariance, x, y)
StationaryRandomFields.spectral_kernel(c::PeriodicCovariance) = PeriodicSpectralKernel(c)

end
