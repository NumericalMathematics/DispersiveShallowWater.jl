# TODO: Only useful as reference for the hyperbolic approximation!
struct SainteMarieEquations1D{Bathymetry <: BathymetryFlat, RealT <: Real} <:
       AbstractSainteMarieEquations{1, 3}
    bathymetry_type::Bathymetry # type of bathymetry
    gravity::RealT # gravitational acceleration
    eta0::RealT # constant still-water surface
end

function SainteMarieEquations1D(; bathymetry_type = bathymetry_flat,
                                  gravity,
                                  eta0 = 0.0)
    SainteMarieEquations1D(bathymetry_type, gravity, eta0)
end

# TODO: more

# TODO: This needs to be implemented.
#       May be found in https://doi.org/10.3934/dcdsb.2015.20.961
# function initial_condition_soliton(x, t, equations::SainteMarieEquations1D, mesh)
#     g = gravity(equations)

#     # setup parameters data
#     h1 = 1.0
#     h2 = 1.2
#     c = sqrt(g * h2)

#     x_t = mod(x - c * t - xmin(mesh), xmax(mesh) - xmin(mesh)) + xmin(mesh)

#     h = h1 + (h2 - h1) * sech(x_t / 2 * sqrt(3 * (h2 - h1) / (h1^2 * h2)))^2
#     v = c * (1 - h1 / h)

#     D = zero(h)
#     b = equations.eta0 - D
#     eta = h + b

#     return SVector(eta, v, D)
# end

