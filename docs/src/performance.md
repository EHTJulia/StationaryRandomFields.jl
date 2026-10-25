# Performance

Native CPU/CUDA and compiled Reactant timings for three synthesis operations on identical resident arrays.
Times are **median ± standard deviation**, in milliseconds per allocating `apply` call.
Standard deviation describes variation across repeated executions, not uncertainty in the median.
GPU speedups in parentheses are ``t_{\mathrm{CPU}}/t_{\mathrm{GPU}}`` using medians for the corresponding implementation: native CPU for CUDA, Reactant CPU for Reactant GPU.

Measured with an AMD EPYC 7543 CPU, NVIDIA A40 GPU, and Julia 1.10.12.
Julia/FFTW use four threads; BLAS uses one. Medians follow warmup and GPU
synchronization, excluding planning, transfers, and compilation. Reactant 0.2.299 uses default compiler and backend threading settings.

## Orthonormal Hartley synthesis

A prepared Hartley transform with orthonormal scaling.

| Precision | Shape | CPU | CUDA | Reactant CPU | Reactant GPU |
| :-- | :-- | --: | --: | --: | --: |
| Float32 | 4096 | 0.0819 ± 0.0366 | 0.0665 ± 0.0299 (1.23×) | 0.0386 ± 0.0215 | 0.0926 ± 0.0430 (0.42×) |
| Float32 | 256 × 256 | 0.5592 ± 0.2168 | 0.0724 ± 0.0422 (7.72×) | 0.7844 ± 0.0928 | 0.0950 ± 0.0870 (8.25×) |
| Float32 | 64 × 64 × 64 | 0.9330 ± 0.2268 | 0.0839 ± 0.0504 (11.13×) | 0.8173 ± 0.3879 | 0.0981 ± 0.0633 (8.33×) |
| Float64 | 4096 | 0.0352 ± 0.0227 | 0.1151 ± 0.0251 (0.31×) | 0.0587 ± 0.0490 | 0.1617 ± 0.0633 (0.36×) |
| Float64 | 256 × 256 | 0.3904 ± 0.1625 | 0.0902 ± 0.0417 (4.33×) | 0.5303 ± 0.2755 | 0.1304 ± 0.0677 (4.07×) |
| Float64 | 64 × 64 × 64 | 3.0929 ± 0.8146 | 0.1636 ± 0.0427 (18.90×) | 1.2555 ± 0.6256 | 0.2077 ± 0.0592 (6.04×) |

## Spectral synthesis

Power-law scaling (index -2, zero DC power) followed by orthonormal Hartley synthesis.

| Precision | Shape | CPU | CUDA | Reactant CPU | Reactant GPU |
| :-- | :-- | --: | --: | --: | --: |
| Float32 | 4096 | 0.0960 ± 0.0271 | 0.0937 ± 0.0339 (1.02×) | 0.0680 ± 0.0381 | 0.0937 ± 0.0444 (0.73×) |
| Float32 | 256 × 256 | 1.6302 ± 0.6049 | 0.1077 ± 0.0450 (15.14×) | 0.7820 ± 0.0741 | 0.1046 ± 0.0890 (7.48×) |
| Float32 | 64 × 64 × 64 | 5.9730 ± 0.6021 | 0.5228 ± 0.0336 (11.42×) | 0.8762 ± 0.3396 | 0.1077 ± 0.0919 (8.14×) |
| Float64 | 4096 | 0.1205 ± 0.0201 | 0.1505 ± 0.0284 (0.80×) | 0.1509 ± 0.0272 | 0.1702 ± 0.0423 (0.89×) |
| Float64 | 256 × 256 | 1.9884 ± 0.5785 | 0.1905 ± 0.0239 (10.44×) | 1.4672 ± 0.1960 | 0.2004 ± 0.0812 (7.32×) |
| Float64 | 64 × 64 × 64 | 11.3058 ± 2.7466 | 0.6933 ± 0.0280 (16.31×) | 2.7813 ± 0.7830 | 0.4848 ± 0.0983 (5.74×) |

## Daubechies 4 synthesis

Periodic orthogonal Db4 with the default maximum depth: 12 levels for 4,096 samples,
8 for 256 × 256, and 6 for 64 × 64 × 64.

| Precision | Shape | CPU | CUDA | Reactant CPU | Reactant GPU |
| :-- | :-- | --: | --: | --: | --: |
| Float32 | 4096 | 0.0962 ± 0.0050 | 0.2304 ± 0.0264 (0.42×) | 0.0466 ± 0.0357 | 0.1610 ± 0.0472 (0.29×) |
| Float32 | 256 × 256 | 2.2958 ± 0.0599 | 0.2566 ± 0.0248 (8.95×) | 1.4601 ± 0.5189 | 0.2017 ± 0.1153 (7.24×) |
| Float32 | 64 × 64 × 64 | 14.1297 ± 0.4340 | 0.4368 ± 0.0285 (32.35×) | 7.6546 ± 1.7138 | 0.2486 ± 0.0806 (30.79×) |
| Float64 | 4096 | 0.0945 ± 0.0091 | 0.2297 ± 0.0298 (0.41×) | 0.0673 ± 0.0268 | 0.1314 ± 0.0532 (0.51×) |
| Float64 | 256 × 256 | 2.5590 ± 0.1028 | 0.2590 ± 0.0231 (9.88×) | 2.6729 ± 0.8136 | 0.2071 ± 0.1092 (12.91×) |
| Float64 | 64 × 64 × 64 | 14.2556 ± 0.4035 | 0.4459 ± 0.0310 (31.97×) | 17.6827 ± 2.8864 | 0.3832 ± 0.0507 (46.15×) |

CUDA uses the Wavelets.jl GPU filter-bank kernels (v0.10.2+); Reactant compiles the separable array implementation.

[Benchmark details and raw results](https://github.com/EHTJulia/StationaryRandomFields.jl/tree/main/benchmark).
