# StationaryRandomFields.jl

StationaryRandomFields represents array-valued latent priors and deterministic maps to
signals. Gaussian latent coefficients can be scaled by a power spectrum and synthesized
with a Hartley transform. Laplace latent coefficients can be synthesized with wavelets.

Start with the [tutorial](tutorial.md), then read the [mathematical conventions](math.md).
The public density interface belongs to the latent prior; signal density evaluation is
outside the current scope.
