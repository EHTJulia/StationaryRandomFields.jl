# Wavelet transforms

A wavelet coefficient controls a feature at a particular location and scale. Synthesis
adds these features to form an image or time series. Laplace coefficients favor a small
number of prominent features compared with an equally scaled Gaussian prior. They are
continuous random variables: this preference does not create exact zeros by itself.
Scale-dependent coefficient weights can favor broad structure over fine texture.

The transform uses periodic orthogonal filter-bank DWT packing compatible
with Wavelets.jl. The approximation occupies the low-index corner; each subsequent
level transforms that corner. By default, the decomposition uses the maximum
number of levels for which every dimension is divisible by ``2^{\mathrm{levels}}``.
An odd-sized axis gives zero levels; an explicit `levels` can request a shallower transform.

For synthesis ``x=SDz``, ``S`` is the inverse wavelet transform and ``D`` is the
diagonal matrix of coefficient weights. A Laplace latent prior
still evaluates `z`. For orthogonal `S` and invertible `D`, this corresponds, up to constants
for fixed weights, to an L1 penalty on ``D^{-1}S^\mathsf{T}x``. Redundant/custom synthesis does not
in general define that analysis penalty on `x`.

See the [Db4 signal](tutorial-wavelet-1d.md) and [image](tutorial-wavelet-2d.md) examples.
