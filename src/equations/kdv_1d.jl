@doc raw"""
    KdVEquation1D(; gravity, D = 1.0, eta0 = 0.0)

Hier noch ganz viel Documentation
"""
struct KdVEquation1D{RealT <: Real} <: AbstractKdVEquation{1, 1}
    gravity::RealT # gravitational acceleration
    D::RealT # still-water depth
    eta0::RealT # constant still-water surface
end

function KdVEquation1D(; gravity, D = 1.0, eta0 = 0.0)
    eta0 == 0.0 || @warn "The still-water surface needs to be 0 for the KdV equations"
    KdVEquation1D(gravity, D, eta0)
end

"""
    initial_condition_convergence_test(x, t, equations::KdVEquation1D, mesh)

A traveling-wave solution used for convergence tests in a periodic domain, here for dimensional variables.
"""
function initial_condition_convergence_test(x, t, equations::KdVEquation1D, mesh)
    g = gravity(equations)
    D = equations.D

    A = 1.0 # amplitude (free parameter)
    K = 1 / 2 * sqrt(3 * A / D^3)
    c = (sqrt(g * D) + A * sqrt(g * D) / (2 * D))

    x_t = mod(x - c * t - xmin(mesh), xmax(mesh) - xmin(mesh)) + xmin(mesh)

    eta = A * sech(K * x_t)^2

    return SVector(eta)
end

# TODO: initial_condition_manufactured
# TODO: source_terms_manufactured

function create_cache(mesh, equations::KdVEquation1D,
                      solver, initial_condition,
                      ::BoundaryConditionPeriodic,
                      RealT, uEltype)
    g = gravity(equations)
    D = equations.D
    c_0 = sqrt(g * D)
    c_1 = sqrt(g / D)

    cache = (; c_0, c_1)
    return cache
end

function rhs!(dq, q, t, mesh, equations::KdVEquation1D, initial_condition,
              ::BoundaryConditionPeriodic, source_terms, solver, cache)
    (; c_0, c_1) = cache
    D = equations.D
    d = D
    g = gravity(equations)
    eta, = q.x
    deta, = dq.x
    deta .= 0.0

    @trixi_timeit timer() "hyperbolic" begin
        deta .= -1.0 * (c_0 * solver.D1 * eta +
                 0.5 * c_1 * eta .* (solver.D1 * eta) +
                 0.5 * c_1 * solver.D1 * (eta .^ 2) +
                 1 / 6 * c_0 * D^2 * (solver.D2 * eta))  #"D2" is the THIRD derivative operator
    end

    #TODO: source terms
    return nothing
end
