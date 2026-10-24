"""
    PeriodicCovariance(grid, model; variance = 1, check_args = true)

Finite periodic covariance with the same paired discrete power as `SpectralScaling`.
`variance` is the variance of each independent latent Gaussian coefficient.
Coordinates are physical positions, with periods `size(grid) .* grid.spacing`.
Off-grid evaluation uses a finite cosine sum; it is not a continuous spectral integral.
Use `prepare(covariance, prototype)` and `apply` for a fast covariance action on the
complete grid. This action does not construct a dense kernel matrix.
"""
struct PeriodicCovariance{G, P, V, F}
    grid::G
    model::P
    variance::V
    frequencies::F
end
function PeriodicCovariance(
        grid::FieldGrid, model::AbstractPowerSpectrumModel; variance = 1, check_args = true
    )
    ndims(grid) == ndims(model) || throw(DimensionMismatch("Spectrum rank mismatch."))
    (!check_args || (isfinite(variance) && variance >= 0)) ||
        throw(ArgumentError("Latent variance must be finite and nonnegative."))
    return PeriodicCovariance(
        grid, model, variance, frequency_axes(grid, typeof(first(grid.spacing)))
    )
end
function PeriodicCovariance(
        grid::FieldGrid, model::AbstractPowerSpectrumModel, latent::UnivariateRandomField{N, D}
    ) where {N, D <: Normal}
    size(grid) == size(latent) || throw(DimensionMismatch("Latent shape mismatch."))
    return PeriodicCovariance(grid, model; variance = var(latent.dist))
end
input_size(c::PeriodicCovariance) = size(c.grid)
output_size(c::PeriodicCovariance) = size(c.grid)

function covariance_term(c, delta, model, index)
    n = ndims(c.grid)
    ν = ntuple(d -> c.frequencies[d][index[d]], n)
    reflected = ntuple(d -> c.frequencies[d][mod(1 - index[d], size(c.grid)[d]) + 1], n)
    return paired_power(model, ν..., reflected...) * cospi(2 * sum(ν .* delta))
end

"""
    covariance_value(covariance, x, y, model = covariance.model)

Evaluate the finite periodic covariance at two coordinate tuples or vectors.
Scalar coordinates are also accepted in one dimension. The optional model argument
allows differentiating updated spectrum parameters without rebuilding grid metadata.
Each evaluation costs O(prod(size(grid))); dense matrix assembly remains a separate
operation from the FHT covariance action.
"""
function covariance_value(c::PeriodicCovariance, x, y, model = c.model)
    n = ndims(c.grid)
    length(x) == length(y) == n || throw(DimensionMismatch("Coordinate rank mismatch."))
    ndims(model) == n || throw(DimensionMismatch("Spectrum rank mismatch."))
    delta = ntuple(d -> x[d] - y[d], n)
    indices = CartesianIndices(size(c.grid))
    total = covariance_term(c, delta, model, first(indices))
    for i in Iterators.drop(indices, 1)
        total += covariance_term(c, delta, model, i)
    end
    return c.variance * total / length(indices)
end
function covariance_value(c::PeriodicCovariance, x::Number, y::Number, model = c.model)
    ndims(c.grid) == 1 || throw(DimensionMismatch("Scalar coordinates require a 1D grid."))
    return covariance_value(c, (x,), (y,), model)
end
(c::PeriodicCovariance)(x, y) = covariance_value(c, x, y)

struct PreparedPeriodicCovariance{H, S, V} <: AbstractPreparedTransform
    hartley::H
    scaling::S
    variance::V
end
input_size(c::PreparedPeriodicCovariance) = input_size(c.hartley)
output_size(c::PreparedPeriodicCovariance) = output_size(c.hartley)
function prepare(c::PeriodicCovariance, prototype)
    return PreparedPeriodicCovariance(
        prepare(HartleyTransform(size(c.grid)), prototype),
        prepare(SpectralScaling(c.grid, c.model), prototype),
        c.variance,
    )
end
apply(c::PreparedPeriodicCovariance, x) = apply(c, x, c.scaling.model)
function apply(c::PreparedPeriodicCovariance, x, model::AbstractPowerSpectrumModel)
    ndims(model) == length(input_size(c)) ||
        throw(DimensionMismatch("Spectrum rank mismatch."))
    coefficients = apply(c.hartley, x)
    power = paired_power.(Ref(model), c.scaling.coordinates..., c.scaling.reflected...)
    weighted = sample_eltype(x).(coefficients .* power .* c.variance)
    return apply(c.hartley, weighted)
end

"""
    spectral_kernel(covariance::PeriodicCovariance)

Wrap a finite periodic covariance as a `KernelFunctions.Kernel`. Load
`KernelFunctions` to activate this optional extension. Generic `kernelmatrix`
evaluation constructs a dense matrix; it does not use the fast grid action or
accelerate GP inference automatically.
"""
function spectral_kernel end
