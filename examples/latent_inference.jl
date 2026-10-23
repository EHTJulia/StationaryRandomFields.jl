using StationaryRandomFields, Distributions, Random

function gaussian_example(rng=MersenneTwister(12))
    dims = (32, 32)
    latent = UnivariateGaussianRandomField(Float64, dims)
    z = rand(rng, latent)
    grid = FieldGrid(dims; spacing=(0.1, 0.1))
    spectrum = SaturatedPowerLaw{2}(-2.0, 0.1, 4.0)
    scaling = prepare(SpectralScaling(grid, spectrum), z)
    synthesis = prepare(HartleyTransform(dims), z)
    signal = apply(synthesis, apply(scaling, z))
    observed = signal .+ 0.1randn(rng, dims)
    function logposterior(z, logpower)
        p = exp(logpower) * spectrum
        prediction = apply(synthesis, apply(scaling, z, p))
        likelihood = -sum(abs2, observed .- prediction) / (2 * 0.1^2)
        return logpdf(latent, z) + logpdf(Normal(), logpower) + likelihood
    end
    return (; latent, z, signal, observed, logposterior)
end

function laplace_example(rng=MersenneTwister(12))
    latent = UnivariateLaplaceRandomField(Float32, (32, 32); b=0.2f0)
    z = rand(rng, latent)
    transform = prepare(WaveletTransform(size(z); levels=3), z)
    signal = apply(transform, z)
    return (; latent, z, signal, logprior=logpdf(latent, z))
end

if abspath(PROGRAM_FILE) == @__FILE__
    example = gaussian_example()
    println("Gaussian latent log posterior: ", example.logposterior(example.z, 0.0))
    println("Laplace latent log prior: ", laplace_example().logprior)
end
