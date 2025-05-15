@doc raw"""
    KdVEquation1D(; gravity, D = 1.0, eta0 = 0.0)

The equation is given by
```math
\begin{aligned}
  \eta_t+\sqrt{g D} \eta_x+3 / 2 \sqrt{g / D} \eta \eta_x+1 / 6 \sqrt{g D} D^2 \eta_{x x x} $= 0.
\end{aligned}
```

and the Semi is:

taken from:

...

Hier noch ganz viel Documentation...
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

"""
    initial_condition_manufactured(x, t, equations::KdVEquation1D, mesh)

A smooth manufactured solution in combination with [`source_terms_manufactured`](@ref).
"""
function initial_condition_manufactured(x, t, equations::KdVEquation1D,
                                        mesh)
    eta = 1 + exp(-t / 2) * sinpi(2 * (x - t / 2))
    return SVector(eta)
end

"""
    source_terms_manufactured(q, x, t, equations::KdVEquation1D, mesh)

A smooth manufactured solution in combination with [`initial_condition_manufactured`](@ref).

calculated using:
```julia
    using Symbolics
    @variables x t 
    Dt = Differential(t)
    Dx = Differential(x)

    @variables g D pi

    η = 1 + exp(-t // 2) * sin(pi*(2 * (x - t // 2)))
    η_t = (Dt(η))
    η_x = (Dx(η))
    η_xxx = (Dx(Dx(Dx(η))))

    source_kdv = η_t + sqrt(g * D) * η_x + 3//2 * sqrt(g / D) * η * η_x + 1//6 * sqrt(g * D) * D^2 * η_xxx

    source_julia_form = simplify(expand_derivatives(source_kdv))
``` 
"""
function source_terms_manufactured(q, x, t, equations::KdVEquation1D)
    g = gravity(equations)
    D = equations.D

    a1 = sinpi(2x - t)
    a2 = cospi(2x - t)
    b1 = exp(-t / 2)
    c0 = sqrt(g * D)
    c1 = sqrt(g / D)

    dq1 = -0.5 * a1 * b1 - pi * a2 * b1 +
          2pi * a2 * c0 * b1 +
          3pi * a2 * c1 * (1 + a1 * b1) * b1 -
          (4 / 3) * D^2 * pi^3 * a2 * c0 * b1

    return SVector(dq1)
end

function create_cache(mesh, equations::KdVEquation1D,
                      solver, initial_condition,
                      ::BoundaryConditionPeriodic,
                      RealT, uEltype)
    g = gravity(equations)
    D = equations.D
    DD = D^2
    c_0 = sqrt(g * D)
    c_1 = 0.5 * sqrt(g / D)

    # We use `DiffCache` from PreallocationTools.jl to enable automatic/algorithmic differentiation
    # via ForwardDiff.jl. 
    # 1: eta
    N = ForwardDiff.pickchunksize(1 * nnodes(mesh))
    template = ones(RealT, nnodes(mesh))

    eta2 = DiffCache(zero(template), N)
    eta2_x = DiffCache(zero(template), N)

    eta_x = DiffCache(zero(template), N)
    eta_xxx = DiffCache(zero(template), N)

    if solver.D1 isa PeriodicUpwindOperators
        D1 = solver.D1.central

        # calculate the third derivative operator using upwind operators
        D3 = sparse(solver.D1.plus) * sparse(solver.D1.minus) * sparse(solver.D1.central)
    else
        D1 = solver.D1
        D3 = solver.D3
    end

    cache = (; D1, D3, c_0, c_1, DD, eta2, eta2_x, eta_x, eta_xxx,)
    return cache
end

function rhs!(dq, q, t, mesh, equations::KdVEquation1D, initial_condition,
              ::BoundaryConditionPeriodic, source_terms, solver, cache)
    eta, = q.x
    deta, = dq.x

    (; D1, D3, c_0, c_1, DD,) = cache
    # In order to use automatic differentiation, we need to extract
    # the storage vectors using `get_tmp` from PreallocationTools.jl
    # so they can also hold dual numbers when needed.
    eta2 = get_tmp(cache.eta2, eta)
    eta2_x = get_tmp(cache.eta2_x, eta2)
    eta_x = get_tmp(cache.eta_x, eta)
    eta_xxx = get_tmp(cache.eta_xxx, eta)

    @trixi_timeit timer() "hyperbolic" begin
        @.. eta2 = eta^2
        # eta2_x = D1 * eta2
        mul!(eta2_x, D1, eta2)

        # eta_x = D1 * eta
        mul!(eta_x, D1, eta)

        @.. deta = -1.0 * (c_0 * eta_x +
                    c_1 * (eta * eta_x +
                           eta2_x)
    end
        
    @trixi_timeit timer() "third-order derivatives" begin
        # eta_xxx = D3 * eta
        mul!(eta_xxx, D3, eta)
        
        @.. deta -= 1 / 6 * c_0 * DD * eta_xxx
    end

    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)

    return nothing
end
