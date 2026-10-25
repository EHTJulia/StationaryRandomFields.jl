# Shared checks run with native CUDA or a compiled Reactant evaluator.
function check_backend_descriptors(to_device, evaluate)
    return @testset "Deferred descriptors on device" begin
        for T in (Float32, Float64), dims in ((16,), (8, 8), (4, 4, 4))
            latent = UnivariateGaussianRandomField(T, dims)
            z = rand(MersenneTwister(17), latent)
            x = to_device(z)
            model = MaternPS{length(dims)}(T(2), T(1.5))
            for spec in (
                    HartleyTransform(), WaveletTransform(wavelet(WT.db4)),
                    TransformChain(SpectralScaling(model), HartleyTransform()),
                    TransformChain(SpectralScaling(model; normalization = :unit_variance), HartleyTransform()),
                    TransformChain(ScaleTransform(T(2)), WaveletTransform(wavelet(WT.db4))),
                )
                p = prepare(spec, latent; prototype = x, spacing = ntuple(_ -> 0.2, length(dims)))
                pc = prepare(spec, latent; spacing = ntuple(_ -> 0.2, length(dims)))
                @test Array(evaluate(p, x)) ≈ apply(pc, z) rtol = T === Float32 ? 1.0e-4 : 1.0e-10
                @test Array(x) == z
                @test !applicable(logpdf, p, x)
            end
            spec = SpectralScaling(model; normalization = :unit_variance)
            p = prepare(spec, latent; prototype = x)
            amplitudes = Array(evaluate(p, to_device(ones(T, dims))))
            @test sum(abs2, amplitudes) / length(amplitudes) ≈ one(T)
        end
    end
end
