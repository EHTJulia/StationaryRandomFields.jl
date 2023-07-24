# Mathematical Foundations

To generate correlated noise, the procedure designated in [Timmer et al. (1995)](https://ui.adsabs.harvard.edu/abs/1995A%26A...300..707T/abstract) is followed. First, random Gaussian noise *G(f)* is generated in the Fourier domain, by drawing from a normal distribution at each point in frequency space. The frequency space grid is computed simply from the size of the signal data itself. Fourier Gaussian noise is then scaled by the square root of a designated Fourier domain power function: 

$$ A({f}) = \sqrt{{P}({f})} = C |{f}|^{_\beta} $$

in which C is a normalization constant and $\beta$ is a negative index. This is the most basic power spectrum available for use, but other special power law functions may also be implemented.

The power-law-scaled Fourier noise, ${F}({f}) = {G}({f}) \cdot {A}({f})$ is Fourier transformed back into the signal domain to produce the output signal noise:

$$
    S(x) = \frac{1}{N} \sum_{f = 0}^{N}F(f)e^{2\pi i f x}.
$$

The above equation is given in one dimension; for data of N dimension, an N-dimensional backward discrete Fourier transform applies.
