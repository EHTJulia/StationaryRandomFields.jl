<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML"></script>

# Mathematical Foundations

To generate correlated noise, the procedure designated in [] is followed. First, random Gaussian noise *G(f)* is generated in the Fourier domain, by drawing from a normal distribution at each point in frequency space. The frequency space grid is computed simply from the size of the signal data itself. Fourier Gaussian noise is then scaled by the square root of a designated Fourier domain power function: 

$$ A({f}) = \sqrt{{P}({f})} = C |{f}|^{_\beta} $$

in which C is a normalization constant and $\beta$ is a negative index. In the case of atmospheric and ISM noise, a saturated power law may be used in which an inner scale $ {f_0} $ and outer scale ${f_1}$ are defined. The saturated power spectrum becomes $C |{f_0}|^{_\beta}$ for all $|{f}|\leq f_0$, and $0$ for $|{f}|\geq f_1$, and remains unchanged in the interval $ f_0 \leq |{f}| \leq f_1 $.


The power-law-scaled Fourier noise, $\bm{F}(\bm{f}) = \bm{G}(\bm{f}) \cdot \bm{A}(\bm{f})$ is Fourier transformed back into the signal domain to produce the output signal noise:

$$
    S(x) = \frac{1}{N} \sum_{f = 0}^{N}F(f)e^{2\pi i f x}.
$$

The above equation is given in one dimension; for data of N dimension, an N-dimensional backward discrete Fourier transform applies.
