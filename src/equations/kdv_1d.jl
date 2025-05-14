@doc raw"""
    KdVEquation1D(; gravity, D = 1.0, eta0 = 0.0)

The equation is given by
```math
\begin{aligned}
  \eta_t + \sqrt{gD}\eta_x + \frac{3}{2}\sqrt{\frac{g}{D}}\eta\eta_x - \frac{1}{6}D^2\eta_{xxt} &= 0.
\end{aligned}
```

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
    DD = D^2
    c_0 = sqrt(g * D)
    c_1 = 0.5 * sqrt(g / D)

    eta2 = zeros(RealT, nnodes(mesh))
    eta2_x = zero(eta2)

    eta_x = zero(eta2)
    eta_xxx = zero(eta2)

    if solver.D1 isa PeriodicUpwindOperators
        D1 = solver.D1.central

        # calculate the third derivative operator using upwind operators
        D3 = sparse(solver.D1.plus) * sparse(solver.D1.minus) * sparse(solver.D1.central)

    else
        D1 = solver.D1
        # The solver struct calls for D2::SecondDerivative, in the case
        # of the KdV equation, one needs a third derivative however
        D3 = solver.D2
    end

    cache = (; D1, D3, eta2, eta2_x, eta_x, eta_xxx, c_0, c_1, DD)
    return cache
end

function rhs!(dq, q, t, mesh, equations::KdVEquation1D, initial_condition,
              ::BoundaryConditionPeriodic, source_terms, solver, cache)
    (; D1, D3, eta2, eta2_x, eta_x, eta_xxx, c_0, c_1, DD) = cache
    D = equations.D
    d = D
    g = gravity(equations)
    eta, = q.x
    deta, = dq.x

    @trixi_timeit timer() "hyperbolic" begin
        @.. eta2 = eta^2
        # eta2_x = D1 * eta2
        mul!(eta2_x, D1, eta2)

        # eta_x = D1 * eta
        mul!(eta_x, D1, eta)

        # eta_xxx = D3 * eta
        mul!(eta_xxx, D3, eta)

        @.. deta = -1.0 * (c_0 * eta_x +
                    c_1 * eta * eta_x +
                    c_1 * eta2_x +
                    1 / 6 * c_0 * DD * eta_xxx)
    end

    #TODO: source terms
    return nothing
end
