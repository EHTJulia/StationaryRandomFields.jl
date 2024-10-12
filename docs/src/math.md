```@meta
CurrentModule = StationaryRandomFields
```

# Mathematical Basis

To generate correlated noise, the procedure designated in [Timmer et al. (1995)](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..707T/abstract) is followed. 
First, random Gaussian noise $G(\vec{f})$ is generated in the Fourier domain, by drawing from a normal distribution at each point in frequency space. 
The frequency space grid is computed simply from the size of the signal data itself. 
Fourier Gaussian noise is then scaled by the square root of a designated Fourier domain power function: 

```math
A(\vec{f}) = \sqrt{ P(\vec{f}) } = \sqrt{ f^{\beta} }
```

in which $\beta$ is a negative index and $f$ denotes the vector norm $|\vec{f}|$. This is the most basic power spectrum available for use, but other special power law functions may also be implemented.

The power-law-scaled Fourier noise, $F(\vec{f}) = G(\vec{f}) \cdot A(\vec{f})$ is inverse Fourier transformed back into the signal domain to produce the output signal noise:

```math
S(\vec{r}) = \frac{1}{N} \sum_{f = 0}^{N} F(\vec{f})e^{2\pi i \vec{f} \cdot \vec{r}}
```

Reversing the process, the approximate power spectrum can be retrieved from a given signal noise $S(\vec{r})$ via a forward Fourier transform:

```math
\bar{P}(\vec{f}) = |F(\vec{f})|^2 = |\frac{1}{N} \sum_{f = 0}^{N} S(\vec{r})e^{-2\pi i \vec{f} \cdot \vec{r}}|^2
```

The approximated power spectrum output by [`get_power_spectrum`](@ref) divides this above result, $\bar{P}(\vec{f})$, by 2 to correct for a positively skewed offset in the squared amplitude of Gaussian noise.  

## Power Laws

StationaryRandomFields.jl implements additional power law functions to be substituted as $P(f)$ in the above equations. [`SaturatedPowerLaw`](@ref) takes the standard form of $P(\vec{f})=f^{\beta}$ only within the interval $f_0 < f < f_1 $, where $f_0$ and $f_1$ are input inner and outer scale cutoffs (`inscale` and `outscale`), respectively. When ${f}\leq f_0$, the spectrum becomes:

```math
P(\vec{f})=f_0^{\beta}
```

and when $f\geq f_1$: 

```math
P(\vec{f})=0.
```


