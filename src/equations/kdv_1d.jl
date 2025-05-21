@doc raw"""
    KdVEquation1D(; gravity, D = 1.0, eta0 = 0.0)

KdV (Korteweg-de Vries) equation in one spatial dimension.
The equation is given by
```math
\begin{aligned}
  \eta_t+\sqrt{g D} \eta_x+3 / 2 \sqrt{g / D} \eta \eta_x+1 / 6 \sqrt{g D} D^2 \eta_{x x x} &= 0.
\end{aligned}
```

The unknown quantity of the KdV equation is the total water height ``\eta``.
The gravitational acceleration `gravity` is denoted by ``g`` and the constant bottom topography (bathymetry) ``b = \eta_0 - D``,
where ``\eta_0`` is the constant still-water surface and ``D`` the still-water depth. The water height above
the bathymetry is therefore given by ``h = \eta - \eta_0 + D``. The KdV equation is only implemented for ``\eta_0 = 0``.

The equations only support a flat bathymetry.

The KdV equation is first introduced by Joseph Valentin Boussinesq (1877) and rediscovered by Diederik Korteweg and Gustav de Vries in 1895.

The semidiscretization implemented here is a modification of the one proposed by Biswas, Ketcheson, Ranocha, and Schütz (2025) for the non-dimensionalized KdV equation ``u_t + u u_x + u_{x x x} = 0.``

The semidiscretization looks the following:
```math
\begin{aligned}
  \eta_t+\sqrt{g D} D_1\eta+ 1 / 2 \sqrt{g / D} \eta D_1 \eta +  1 / 2 \sqrt{g / D} D_1 \eta^2 +1 / 6 \sqrt{g D} D^2 D_3\eta &= 0.
\end{aligned}
```
where ``D_1`` is a first-derivative operator, ``D_3`` a third-derivative operator, and ``D`` the still-water depth.

It conserves
- the total water mass (integral of ``\eta``) as a linear invariant
and if upwind operators (``D_3 = D_{1,+} D_1 D_{1,-}``) or wide-stencil operators (``D_3 = D_1^3``) are used for the third derivative, it also conserves
- the energy (integral of ``1/2\eta^2``)

for periodic boundary conditions.

- Diederik Korteweg and Gustav de Vries (1895)
  On the change of form of long waves advancing in a rectangular canal, and on a new type of long stationary waves
  [DOI: 10.1080/14786449508620739](https://doi.org/10.1080/14786449508620739)
- Abhijit Biswas, David I. Ketcheson, Hendrik Ranocha and Jochen Schütz (2025)
  Traveling-Wave Solutions and Structure-Preserving Numerical Methods for a Hyperbolic Approximation of the Korteweg-de Vries Equation
  [DOI: 10.1007/s10915-025-02898-x](https://doi.org/10.1007/s10915-025-02898-x)
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
    c0 = sqrt(g * D)
    c = 1.5 * c0
    A = 2 * D * (c - c0) / c0
    K = 1 / 2 * sqrt(3 * A / D^3)
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

How it was calculated, is described in:
https://github.com/NumericalMathematics/DispersiveShallowWater.jl/pull/198#discussion_r2090805751
"""
function source_terms_manufactured(q, x, t, equations::KdVEquation1D)
    g = gravity(equations)
    D = equations.D

    a1 = sinpi(2x - t)
    a2 = cospi(2x - t)
    b1 = exp(-t / 2)
    c0 = sqrt(g * D)
    c1 = sqrt(g / D)

    s1 = -0.5 * a1 * b1 - pi * a2 * b1 +
         2pi * a2 * c0 * b1 +
         3pi * a2 * c1 * (1 + a1 * b1) * b1 -
         (4 / 3) * D^2 * pi^3 * a2 * c0 * b1

    return SVector(s1)
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
    # nvariables(equations) = 1: eta
    N = ForwardDiff.pickchunksize(nvariables(equations) * nnodes(mesh))
    template = ones(RealT, nnodes(mesh))

    tmp_1 = DiffCache(zero(template), N)
    tmp_2 = DiffCache(zero(template), N)

    cache = (; c_0, c_1, DD, tmp_1, tmp_2)
    return cache
end

function rhs!(dq, q, t, mesh, equations::KdVEquation1D, initial_condition,
              ::BoundaryConditionPeriodic, source_terms, solver, cache)
    eta, = q.x
    deta, = dq.x

    (; c_0, c_1, DD) = cache
    # In order to use automatic differentiation, we need to extract
    # the storage vectors using `get_tmp` from PreallocationTools.jl
    # so they can also hold dual numbers when needed.
    tmp_1 = get_tmp(cache.tmp_1, eta)
    tmp_2 = get_tmp(cache.tmp_2, eta)

    @trixi_timeit timer() "third-order derivatives" begin
        if solver.D1 isa PeriodicUpwindOperators && isnothing(solver.D3)
            # eta_xxx = Dp * Dc * Dm * eta
            mul!(tmp_1, solver.D1.minus, eta)
            mul!(tmp_2, solver.D1.central, tmp_1)
            mul!(tmp_1, solver.D1.plus, tmp_2)

            # set D1 for hyperbolic terms
            D1 = solver.D1.central
        else
            # eta_xxx = D3 * eta
            mul!(tmp_1, solver.D3, eta)

            # set D1 for hyperbolic terms
            D1 = solver.D1
        end

        # deta = 1 / 6 sqrt(g * D) D^2 eta_xxx
        @.. deta = -1 / 6 * c_0 * DD * tmp_1
    end

    @trixi_timeit timer() "hyperbolic" begin
        # eta2 = eta^2
        @.. tmp_1 = eta^2
        # eta2_x = D1 * eta2
        mul!(tmp_2, D1, tmp_1)

        # eta_x = D1 * eta
        mul!(tmp_1, D1, eta)

        # deta -= sqrt(g * D) * eta_x + 1 / 2 * sqrt(g / D) * (eta * eta_x + eta2_x)
        @.. deta -= (c_0 * tmp_1 + c_1 * (eta * tmp_1 + tmp_2))
    end

    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)

    return nothing
end
