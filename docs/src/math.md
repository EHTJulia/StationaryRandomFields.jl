# Mathematical Basis

## The idea before the equations

Choose how much variation you expect on each scale, draw independent random
coefficients, and turn those coefficients into a signal. For Gaussian synthesis,
large-scale modes get larger amplitudes when the power spectrum is steeper. The
[illustrated tutorial](tutorial.md) shows each step before introducing any inference.

This follows the spectral-synthesis idea of **Timmer & König (1995),
“On Generating Power Law Noise,” A&A 300, 707–710**
([SciX](https://scixplorer.org/abs/1995A%26A...300..707T/abstract)). Their algorithm
randomizes both Fourier amplitudes and phases. Here independent Gaussian Hartley
coefficients give a real-valued representation of the same Gaussian spectral-synthesis
principle on a finite periodic grid, with the multidimensional normalization below.
We do not fix every realized Fourier amplitude to the square root of its target power:
the power spectrum is an ensemble expectation, and individual realizations fluctuate.

In the equations, ``z`` means the latent coefficient array, ``x`` the synthesized signal,
``M`` the number of array elements, and ``P`` the power assigned to a discrete mode. A
covariance describes how values at two signal locations vary together across draws.

## Choose a topic

- [Latent fields and inference](math-latent.md): distributions, sampling, and prior densities.
- [Spectral transforms](math-spectral.md): Gaussian synthesis, discrete power, and covariance.
- [Wavelet transforms](math-wavelet.md): multiscale synthesis and Laplace coefficient priors.
