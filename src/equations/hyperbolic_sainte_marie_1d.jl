@doc raw"""
    HyperbolicSainteMarieEquations1D(; bathymetry_type = bathymetry_mild_slope,
                                       gravity,
                                       eta0 = 0.0,
                                       alpha = 3.0)

Hyperbolic approximation of the Sainte-Marie system
[`SainteMarieEquations1D`](@ref) in one spatial dimension derived by
Escalante, Dumbser and Castro (2019).
The equations are given by
```math
\begin{aligned}
  h_t + (h v)_x &= 0,\\
  h v_t + \frac{1}{2} g (h^2)_x + \frac{1}{2} h (v^2)_x
    + (h p)_x &= -(g h + 2 p) b_x,\\
  h w_t + h v w_x &= 2 p,\\
  h p_t + h v p_x + h c^2 \bigl(v_x + (w - v b_x) / (h / 2)\bigr) &= 0.
\end{aligned}
```
The unknown quantities of the Sainte-Marie equations are the
total water height ``\eta = h + b`` and the velocity ``v``.
The gravitational acceleration `gravity` is denoted by ``g`` and the bottom topography
(bathymetry) ``b = \eta_0 - D``. The water height above the bathymetry
is therefore given by ``h = \eta - \eta_0 + D``.
The total water height is therefore given by ``\eta = h + b``.

There are two additional variables compared to the [`SainteMarieEquations1D`](@ref):
the auxiliary variable ``w \approx -h v_x / 2 + v b_x`` and the non-hydrostatic pressure ``p``.
In the formal limit ``c \to \infty``, the hyperbolic approximation recovers the original Sainte-Marie system.
The hyperbolization parameter is determined as ``c = \alpha \sqrt{g \eta_0}`` set by the keyword arguments `alpha` (``\alpha``), `gravity` (``g``), and `eta0` (``\eta_0``).
The larger the value of ``\alpha``, the better the approximation of the original system, but also the stiffer the system.

!!! note "Initial conditions"
    The `HyperbolicSainteMarieEquations1D` allow two options for
    specifying initial conditions:
    1. Returning the full set of variables `q = (η, v, D, w, p)`
    2. Returning a reduced set of variables `q = (η, v, D)` as required
       for the limit system [`SainteMarieEquations1D`](@ref) of the
       hyperbolic approximation. The remaining variables `w` and `p` are
       initialized using the default initialization ``w \approx -h v_x / 2 + v b_x``
       and ``p \approx 0`` using the derivative operator of the `solver`.

Two types of `bathymetry_type` are supported:
- [`bathymetry_flat`](@ref): flat bathymetry (typically ``b = 0`` everywhere)
- [`bathymetry_mild_slope`](@ref): variable bathymetry with mild-slope approximation

References for the Sainte-Marie system and its hyperbolization can be found in
- Sainte-Marie (2011)
  Vertically averaged models for the free surface non-hydrostatic
  Euler system: derivation and kinetic interpretation
  [DOI: 10.1142/S0218202511005118](https://doi.org/10.1142/S0218202511005118)
- Bristeau, Mangeney, Sainte-Marie and Seguin (2015)
  An energy-consistent depth-averaged Euler system:
  derivation and properties
  [DOI: 10.3934/dcdsb.2015.20.961](https://doi.org/10.3934/dcdsb.2015.20.961)
- Aïssiouene, Bristeau, Godlewski, Mangeney, Parés Madroñal and Sainte-Marie (2020)
  A two-dimensional method for a family of dispersive shallow water models
  [DOI: 10.5802/smai-jcm.66](https://doi.org/10.5802/smai-jcm.66)
- Escalante, Dumbser and Castro (2019)
  An efficient hyperbolic relaxation system for dispersive non-hydrostatic
  water waves and its solution with high order discontinuous Galerkin schemes
  [DOI: 10.1016/j.jcp.2019.05.035](https://doi.org/10.1016/j.jcp.2019.05.035)

The energy-conserving semidiscretization is basically taken from the
following reference, but adapted to the standard style of
DispersiveShallowWater.jl to use the primitive variables `q = (η, v, D, w, p)`.
- Artiano and Ranocha (2026)
  On Affordable High-Order Entropy-Conservative/Stable and
  Well-Balanced Methods for Nonconservative Hyperbolic Systems
  [DOI: 10.48550/arXiv.2603.18978](https://arxiv.org/abs/2603.18978)
"""
struct HyperbolicSainteMarieEquations1D{Bathymetry <:
                                        Union{BathymetryFlat, BathymetryMildSlope},
                                        RealT <: Real} <: AbstractSainteMarieEquations{1, 5}
    bathymetry_type::Bathymetry # type of bathymetry
    gravity::RealT # gravitational acceleration
    eta0::RealT # constant still-water surface
    alpha::RealT # hyperbolic relaxation parameter (→ ∞ for Sainte-Marie)
    c_squared::RealT # c^2 = alpha^2 * g * eta0
end

function HyperbolicSainteMarieEquations1D(; bathymetry_type = bathymetry_mild_slope,
                                          gravity,
                                          eta0,
                                          alpha = 3.0)
    c_squared = alpha^2 * gravity * eta0
    return HyperbolicSainteMarieEquations1D(bathymetry_type, gravity, eta0, alpha,
                                            c_squared)
end

# TODO: It would also make sense to use the conservative form as default.
function varnames(::typeof(prim2prim), ::HyperbolicSainteMarieEquations1D)
    return ("η", "v", "D", "w", "p")
end
function varnames(::typeof(prim2cons), ::HyperbolicSainteMarieEquations1D)
    return ("h", "hv", "b", "hw", "hp")
end

"""
    prim2phys(q, equations::HyperbolicSainteMarieEquations1D)

Return the physical variables ``\\eta, v, D`` used also by the
[`SainteMarieEquations1D`](@ref) from the main variables
`q` for the [`HyperbolicSainteMarieEquations1D`](@ref).
"""
function prim2phys(q, ::HyperbolicSainteMarieEquations1D)
    eta, v, D = q
    return SVector(eta, v, D)
end

varnames(::typeof(prim2phys), ::HyperbolicSainteMarieEquations1D) = ("η", "v", "D")

is_hyperbolic_appproximation(::HyperbolicSainteMarieEquations1D) = Val{true}()

function hyperbolic_approximation_limit(equations::HyperbolicSainteMarieEquations1D)
    (; bathymetry_type, gravity, eta0) = equations
    return SainteMarieEquations1D(bathymetry_type, gravity, eta0)
end

# compute the auxiliary variables from the main ones
# using the default initialization
function set_approximation_variables!(q, mesh,
                                      equations::HyperbolicSainteMarieEquations1D,
                                      solver)
    (; D1) = solver
    eta, v, D, w, p = q.x

    # This is not true but just a simple initialization used also in
    # at least one experiment of Escalante, Dumbser, Castro (2019).
    # In general, the non-hydrostatic pressure involves also time derivatives
    # of the velocity.
    @.. p = 0

    # w ≈ -(1/2) h v_x + v b_x
    # with h = (h + b) + (eta0 - b) - eta0
    mul!(w, D1, v)
    @.. w = -0.5 * (eta + D - equations.eta0) * w
    if !(equations.bathymetry_type isa BathymetryFlat)
        b = equations.eta0 .- D
        b_x = D1 * b
        @.. w = w + v * b_x
    end

    return nothing
end

"""
    initial_condition_manufactured(x, t, equations::HyperbolicSainteMarieEquations1D, mesh)

A smooth manufactured solution in combination with
[`source_terms_manufactured`](@ref).
"""
function initial_condition_manufactured(x, t,
                                        equations::HyperbolicSainteMarieEquations1D,
                                        mesh)
    eta = 2 + cospi(2 * (x - 2 * t))
    b = -5 - 2 * cospi(2 * x)
    # h = eta - b = 7 + 2 * cospi(2 * x) + cospi(2 * (x - 2 * t))
    v = sinpi(2 * (x - t / 2))
    D = equations.eta0 - b
    # w = -h v_x / 2 + v b_x
    w = -pi * cospi(t - 2 * x) * (7 + 2 * cospi(2 * x) + cospi(2 * (x - 2 * t))) -
        4 * pi * sinpi(t - 2 * x) * sinpi(2 * x)
    p = cospi(2 * x - 3 * t)
    return SVector(eta, v, D, w, p)
end

"""
    source_terms_manufactured(q, x, t, equations::HyperbolicSainteMarieEquations1D)

A smooth manufactured solution in combination with
[`initial_condition_manufactured`](@ref).
"""
function source_terms_manufactured(q, x, t,
                                   equations::HyperbolicSainteMarieEquations1D)
    g = gravity(equations)
    pi_2 = pi^2

    a1 = cospi(t - 4 * x)
    a2 = cospi(t - 2 * x)
    a3 = cospi(5 * t - 4 * x)
    a4 = sinpi(4 * t - 2 * x)
    a5 = cospi(3 * t - 2 * x)
    a6 = sinpi(2 * x)
    a7 = cospi(2 * x)
    a8 = sinpi(t - 2 * x)
    a9 = sinpi(3 * t - 2 * x)
    a10 = cospi(4 * t - 2 * x)
    a11 = sinpi(4 * x - 2 * t)
    a12 = sinpi(t - 4 * x)
    a13 = sinpi(5 * t - 4 * x)
    a14 = (7 + 2 * a7 + a10)

    # Source terms for variable bathymetry
    s1 = 2 * pi * (2 * a1 + 7 * a2 + a3 - 2 * a4)
    s2 = pi * (-a2 + a5 * (4 * a6 + 2 * a4) / a14 + 2 * g * a4 + a11 + 2 * a9)
    s3 = zero(s1)
    s4 = -2 * a5 / a14 +
         pi_2 * (14 * a8^2 - 4 * a2 * (a6 - a4) + a8 * (a14 + 12 * a12 + 2 * a13))
    s5 = -pi * (3 + 2 * a8) * a9

    return SVector(s1, s2, s3, s4, s5)
end

dingemans_calibration(equations::HyperbolicSainteMarieEquations1D) = 1.0

function create_cache(mesh, equations::HyperbolicSainteMarieEquations1D,
                      solver, initial_condition,
                      boundary_conditions::BoundaryConditionPeriodic,
                      RealT, uEltype)
    # We use `DiffCache` from PreallocationTools.jl to enable automatic/algorithmic differentiation
    # via ForwardDiff.jl. We also pass the second argument determining the chunk size since the
    # typical use case is to compute Jacobians of the full `rhs!` evaluation, where the complete
    # state vector is `q`, which is bigger than the storage for a single scalar variable.
    # nvariables(equations) = 5: eta, v, D, w, H
    N = ForwardDiff.pickchunksize(nvariables(equations) * nnodes(mesh))
    template = ones(RealT, nnodes(mesh))
    h = DiffCache(template, N)
    b = DiffCache(zero(template), N)
    b_x = DiffCache(zero(template), N)
    h_x = DiffCache(zero(template), N)
    v_x = DiffCache(zero(template), N)
    hv_x = DiffCache(zero(template), N)
    v2_x = DiffCache(zero(template), N)
    h_hpb_x = DiffCache(zero(template), N)
    w_x = DiffCache(zero(template), N)
    hvw_x = DiffCache(zero(template), N)
    p_x = DiffCache(zero(template), N)
    hp_x = DiffCache(zero(template), N)
    hvp_x = DiffCache(zero(template), N)
    tmp = DiffCache(zero(template), N)

    cache = (; h, b, b_x, h_x, v_x, hv_x, v2_x, h_hpb_x,
             w_x, hvw_x, p_x, hp_x, hvp_x, tmp)
    return cache
end

# Discretization that conserves
# - the total water mass (integral of ``h``) as a linear invariant
# - the total momentum/discharge as a quadratic invariant for periodic BCs
# - the total modified energy
function rhs!(dq, q, t, mesh,
              equations::HyperbolicSainteMarieEquations1D,
              initial_condition,
              boundary_conditions::BoundaryConditionPeriodic,
              source_terms,
              solver, cache)
    # Unpack physical parameters and SBP operator `D1`
    g = gravity(equations)
    (; c_squared, bathymetry_type) = equations
    (; D1) = solver

    # `q` and `dq` are `ArrayPartition`s. They collect the individual
    # arrays for the total water height `eta = h + b`, the velocity `v`,
    # and the additional variables `w` and `p`.
    eta, v, D, w, p = q.x
    dh, dv, dD, dw, dp = dq.x # dh = deta since b is constant in time
    fill!(dD, zero(eltype(dD)))

    @trixi_timeit timer() "hyperbolic terms" begin
        # Compute all derivatives required below

        # First, we extract temporary storage from the `cache`.
        # Since we use `DiffCache` from PreallocationTools.jl, we need to extract the
        # appropriate arrays using `get_tmp` and need to pass an array with the element
        # type we want to use, e.g., plain `Float64` or some dual numbers when using AD.
        h = get_tmp(cache.h, eta)
        b = get_tmp(cache.b, eta)
        b_x = get_tmp(cache.b_x, eta)
        h_x = get_tmp(cache.h_x, eta)
        v_x = get_tmp(cache.v_x, eta)
        hv_x = get_tmp(cache.hv_x, eta)
        v2_x = get_tmp(cache.v2_x, eta)
        h_hpb_x = get_tmp(cache.h_hpb_x, eta)
        w_x = get_tmp(cache.w_x, eta)
        hvw_x = get_tmp(cache.hvw_x, eta)
        p_x = get_tmp(cache.p_x, eta)
        hp_x = get_tmp(cache.hp_x, eta)
        hvp_x = get_tmp(cache.hvp_x, eta)
        tmp = get_tmp(cache.tmp, eta)

        @.. b = equations.eta0 - D
        @.. h = eta - b
        if !(bathymetry_type isa BathymetryFlat)
            mul!(b_x, D1, b)
        end

        # h_x = D1 * h
        mul!(h_x, D1, h)

        # v_x = D1 * v
        mul!(v_x, D1, v)

        # hv_x = D1 * (h * v)
        @.. tmp = h * v
        mul!(hv_x, D1, tmp)

        # v2_x = D1 * (v.^2)
        @.. tmp = v^2
        mul!(v2_x, D1, tmp)

        # h_hpb_x = D1 * (h .* eta)
        @.. tmp = h * eta
        mul!(h_hpb_x, D1, tmp)

        # w_x = D1 * w
        mul!(w_x, D1, w)

        # hvw_x = D1 * (h * v * w)
        @.. tmp = h * v * w
        mul!(hvw_x, D1, tmp)

        # p_x = D1 * p
        mul!(p_x, D1, p)

        # hp_x = D1 * (h * p)
        @.. tmp = h * p
        mul!(hp_x, D1, tmp)

        # hvp_x = D1 * (h * v * p)
        @.. tmp = h * v * p
        mul!(hvp_x, D1, tmp)

        # Plain: h_t + (h v)_x = 0
        #
        # Split form for energy conservation:
        # h_t + h_x v + h v_x = 0
        @.. dh = -(h_x * v + h * v_x)
        if boundary_conditions isa BoundaryConditionReflecting
            # dh[1] -= h[1] * v[1] / left_boundary_weight(D1)
            # dh[end] += h[end] * v[end] / right_boundary_weight(D1)
            # FIXME
            error("Reflecting boundary conditions are not implemented yet.")
        end

        # Plain: h v_t + h v v_x + g (h + b) h_x
        #              + ... = 0
        #
        # Split form for energy conservation:
        # h v_t + g (h (h + b))_x - g (h + b) h_x
        #       + 1/2 h (v^2)_x - 1/2 v^2 h_x  + 1/2 v (h v)_x - 1/2 h v v_x
        #       + (h p)_x + 2 p b_x = 0
        @.. dv = -(g * h_hpb_x - g * eta * h_x
                   + 0.5 * (h * v2_x - v^2 * h_x)
                   + 0.5 * v * (hv_x - h * v_x)
                   + hp_x + 2 * p * b_x) / h

        # Plain: h w_t + h v w_x = 2 p
        #
        # Split form for energy conservation:
        # h w_t + 1/2 (h v w)_x + 1/2 h v w_x
        #       - 1/2 h_x v w - 1/2 h w v_x = 2 p
        @.. dw = (-0.5 * (hvw_x + h * v * w_x + dh * w) + 2 * p) / h

        # Plain: h p_t + ...
        #
        # Split form for energy conservation:
        # h p_t + 1/2 (h v p)_x + 1/2 h v p_x
        #       - 1/2 h_x v p - 1/2 h v_x p
        #       + c^2 h v_x + 2 c^2 w
        #       - 2 c^2 v b_x = 0
        @.. dp = (-0.5 * (hvp_x + h * v * p_x + dh * p)
                  -
                  c_squared * h * v_x - 2 * c_squared * w
                  +
                  2 * c_squared * v * b_x) / h
    end

    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)

    return nothing
end

@inline function prim2cons(q, equations::HyperbolicSainteMarieEquations1D)
    h = waterheight(q, equations)
    v = velocity(q, equations)
    b = bathymetry(q, equations)

    hv = h * v
    hw = h * q[4]
    hp = h * q[5]
    return SVector(h, hv, b, hw, hp)
end

@inline function cons2prim(u, equations::HyperbolicSainteMarieEquations1D)
    h, hv, b, hw, hp = u

    eta = h + b
    v = hv / h
    D = equations.eta0 - b
    w = hw / h
    p = hp / h
    return SVector(eta, v, D, w, p)
end

# The entropy/energy takes the whole `q` for every point in space
function energy_total_modified!(e, q_global,
                                equations::HyperbolicSainteMarieEquations1D,
                                cache)
    # unpack physical parameters and SBP operator `D1`
    g = gravity(equations)
    (; c_squared) = equations
    h = get_tmp(cache.h, q_global)
    b = get_tmp(cache.b, q_global)

    # `q_global` is an `ArrayPartition`. It collects the individual arrays for
    # the total water height `eta = h + b` and the velocity `v`.
    eta, v, D, w, p = q_global.x
    @.. b = equations.eta0 - D
    @.. h = eta - b

    # 1/2 g eta^2 + 1/2 h v^2 + 1/2 h w^2 + 1/(2 c^2) h p^2
    @.. e = 1 / 2 * g * eta^2 + 1 / 2 * h * v^2 + 1 / 2 * h * w^2 +
            1 / (2 * c_squared) * h * p^2

    return e
end
