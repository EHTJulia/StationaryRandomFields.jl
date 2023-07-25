# Mathematical Basis

To generate correlated noise, the procedure designated in [Timmer et al. (1995)](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..707T/abstract) is followed. First, random Gaussian noise *G(f)* is generated in the Fourier domain, by drawing from a normal distribution at each point in frequency space. The frequency space grid is computed simply from the size of the signal data itself. Fourier Gaussian noise is then scaled by the square root of a designated Fourier domain power function: 

$$ A({f}) = \sqrt{{P}({f})} = \sqrt{C |{f}|^{\beta}} $$

in which C is a normalization constant and $\beta$ is a negative index. This is the most basic power spectrum available for use, but other special power law functions may also be implemented.

The power-law-scaled Fourier noise, ${F}({f}) = {G}({f}) \cdot {A}({f})$ is Fourier transformed back into the signal domain to produce the output signal noise:

$$
    S(x) = \frac{1}{N} \sum_{f = 0}^{N}F(f)e^{2\pi i f x}.
$$

The above equation is given in one dimension; for data of N dimension, an N-dimensional backward discrete Fourier transform applies.

## Power Laws

StationaryRandomFields.jl implements additional power law functions to be substituted as $P(f)$ in the above equations. `SaturatedPowerLaw` takes the standard form of $P(f)=|{f}|^{\beta}$ only within the interval $f_0 \leq |{f}| \leq f_1 $, where $f_0$ and $f_1$ are input inner and outer scales, respectively. The spectrum becomes $C |{f_0}|^{\beta}$ for all $|{f}|\leq f_0$, and $0$ for $|{f}|\geq f_1$.

`PhaseScreenPowerLaw` is specific to Interstellar-Medium-generated fluctuations, intended for use in Event Horizon Telescope imaging of Sagittarius A*. It describes the 2-dimensional phase screen approximation of the ISM, following the equation:

$$ P(\vec{f} ) = \bar{Q} \cdot (|\vec{f}| r_{\text{in}})^{-(\alpha + 2)} \cdot e^{-(|\vec{f}| r_{\text{in}})^2} \cdot P_{\phi}(\phi)$$

in which $\phi$ is the angular polar coordinate of $\vec{f}$. The constant $\bar{Q}$, inner scale $r_{in}$, index $\alpha$, and function $P_{\phi}(\phi)$ are taken as input parameters. This type is intended for usage in tandem with the (under-development) package `ScatteringOptics.jl`.

The amplitude-mapping function for `PhaseScreenPowerLaw` takes an optional parameter for time in hours, `t_hr`. The function shifts in Fourier space to simulate a phase screen moving according to the screen's velocity by scaling amplitude $A(f)$ with the term:
$$\text{exp}({\frac{2π i}{N} (s v_x \tau + t v_y \tau)})$$

in which s and t represent the cartesian coordinates of frequency $\vec{f}$ in Fourier space and $\tau$ is the given time. $v_x$ and $v_y$ are the x and y components of the phase screen velocity in physical space, optionally defined in the construction of the `PhaseScreenPowerLaw` object as parameters `Vx_km_per_s` and `Vy_km_per_s`.
