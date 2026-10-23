# Introduction

The model separates three responsibilities:

 1. `UnivariateRandomField` stores the latent shape and a shared continuous scalar distribution.
 2. An `AbstractFieldTransform` describes a deterministic computation and its shapes.
 3. `prepare(transform, prototype)` creates reusable backend-specific execution metadata.

`TransformedRandomField` composes a latent distribution and a transform for sampling.
Its output shape may differ from the latent shape. Custom transforms implement
`input_size`, `output_size`, `prepare`, and `apply`; `apply!` has a copying fallback.
Built-in Hartley, spectral, identity, and wavelet transforms preserve shape.

Array storage selects the backend. Coefficient-weight arrays must be on the same backend as their input; transfer trainable weights explicitly.
Distribution objects do not contain a GPU device or
an RNG. `rand(field)` follows the Distributions CPU-allocation convention; for device-prepared fields, use `rand!(rng, field, device_output)`. CPU sampling accepts an explicit `AbstractRNG`. CUDA sampling should use a CUDA
RNG or `rand!(field, cuarray)`. Compiled Reactant sampling should use an explicit
`ReactantRNG`; keep RNG state outside differentiation.

Preparation can allocate, copy static frequency metadata to a device, and build plans.
Prepare outside a differentiated or compiled objective. Wavelet filters are copied at preparation. Model-dependent spectral
amplitudes are recomputed during execution. Prepared objects may be reused with matching
array types and shapes. Prepare the original descriptor separately for each concurrent task and use separate output buffers. Hartley
`apply!` rejects overlapping input/output storage.

The package name is retained. Stationarity proofs, shift-invariance tests, and
KernelFunctions integration are outside this redesign's scope.
