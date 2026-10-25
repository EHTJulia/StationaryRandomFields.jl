# Display model parameters without exposing storage type parameters.
function Base.show(io::IO, p::Union{ConstantPowerSpectrum, SinglePowerLaw, SaturatedPowerLaw, MaternPS, MarkovPS})
    print(io, nameof(typeof(p)), "{", ndims(p), "}(")
    for (i, key) in enumerate(fieldnames(typeof(p)))
        i > 1 && print(io, ", ")
        print(io, key, "=")
        show(io, getfield(p, key))
    end
    return print(io, ")")
end
function Base.show(io::IO, p::SumPowerSpectrum)
    print(io, "SumPowerSpectrum(")
    join(io, p.models, ", ")
    return print(io, ")")
end
function Base.show(io::IO, p::ModifiedPowerSpectrumModel)
    print(io, "modify_ps(")
    show(io, p.model)
    for modifier in p.modifiers
        print(io, ", ")
        show(io, modifier)
    end
    return print(io, ")")
end
function Base.show(io::IO, m::PowerSpectrumModifier)
    print(io, nameof(typeof(m)), "(")
    for (i, key) in enumerate(fieldnames(typeof(m)))
        i > 1 && print(io, ", ")
        print(io, key, "=")
        show(io, getfield(m, key))
    end
    return print(io, ")")
end
function Base.show(io::IO, t::Union{HartleyTransform, PreparedHartley})
    print(io, "HartleyTransform(size=", t.dims, ", normalization=orthonormal")
    t isa PreparedHartley && print(io, ", prepared=true")
    return print(io, ")")
end
function Base.show(io::IO, t::Union{SpectralScaling, PreparedSpectralScaling})
    print(io, "SpectralScaling(size=", input_size(t), ", model=")
    show(io, t.model)
    print(io, ", action=√P scaling")
    t isa SpectralScaling && print(io, ", spacing=", t.grid.spacing)
    t isa PreparedSpectralScaling && print(io, ", prepared=true")
    return print(io, ")")
end
function Base.show(io::IO, t::Union{WaveletTransform, PreparedWavelet})
    print(
        io, "WaveletTransform(size=", t.dims, ", basis=", Wavelets.WT.name(t.wavelet),
        ", levels=", t.levels, ", boundary=periodic, direction=synthesis"
    )
    t isa PreparedWavelet && print(io, ", prepared=true")
    return print(io, ")")
end
function Base.show(io::IO, t::IdentityTransform)
    return print(io, "IdentityTransform(size=", t.dims, ")")
end
function Base.show(io::IO, t::ScaleTransform)
    print(io, "ScaleTransform(size=", t.dims, ", scale=")
    t.scale isa AbstractArray ? summary(io, t.scale) : show(io, t.scale)
    return print(io, ")")
end
function Base.show(io::IO, t::TransformChain)
    print(io, "TransformChain(")
    for (i, part) in enumerate(t.transforms)
        i > 1 && print(io, ", ")
        show(io, part)
    end
    return print(io, ")")
end
function Base.show(io::IO, ::MIME"text/plain", t::TransformChain)
    print(io, "TransformChain")
    has_shape(t) && print(io, " (", join(input_size(t), " × "), ")")
    for (i, part) in enumerate(t.transforms)
        print(io, "\n  ", i, ". ")
        show(io, part)
    end
    return
end
