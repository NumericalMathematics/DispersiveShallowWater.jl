@doc raw"""
    HyperbolicSerreGreenNaghdiEquations1D(; bathymetry_type = bathymetry_mild_slope,
                                          gravity,
                                          eta0 = 0.0,
                                          lambda)

Hyperbolic approximation of the Serre-Green-Naghdi system in one spatial
dimension. The equations for flat bathymetry are given by
```math
\begin{aligned}
  h_t + (h v)_x &= 0,\\
  h v_t + \frac{1}{2} g (h^2)_x + \frac{1}{2} h (v^2)_x
    + \biggl( \frac{\lambda}{3} H (1 - H / h) \biggr)_x &= 0,\\
  h w_t + h v w_x &= \lambda (1 - H / h),\\
  H_t + H_x u &= w.
\end{aligned}
```
The unknown quantities of the hyperbolized Serre-Green-Naghdi equations are the
total water height ``\eta = h + b`` and the velocity ``v``.
The gravitational acceleration `gravity` is denoted by ``g`` and the bottom topography
(bathymetry) ``b = \eta_0 - D``. The water height above the bathymetry
is therefore given by ``h = \eta - \eta_0 + D``.
The total water height is therefore given by ``\eta = h + b``.

There are two additional variables ``w \approx -h v_x`` and ``H \approx h``
compared to the [`SerreGreenNaghdiEquations1D`](@ref). In the original papers
of Gavrilyuk et al., the variable ``H`` is called ``\eta``. Here, we
use ``\eta`` for the total water height and ``H`` for auxiliary variable
introduced in the hyperbolic approximation.

!!! note "Initial conditions"
    The `HyperbolicSerreGreenNaghdiEquations1D` allow two options for
    specifying initial conditions:
    1. Returning the full set of variables `q = (η, v, D, w, H)`
    2. Returning a reduced set of variables `q = (η, v, D)` as required
       for the limit system [`SerreGreenNaghdiEquations1D`](@ref) of the
       hyperbolic approximation. The remaining variables `w` and `H` are
       initialized using the default initialization ``w \approx -h v_x``
       and ``H \approx h`` using the derivative operator of the `solver`.

The relaxation parameter `lambda` (``\lambda``) introduced to obtain
this hyperbolic approximation of the [`SerreGreenNaghdiEquations1D`](@ref)
influences the stiffness of the system. For ``\lambda \to \infty``, the
hyperbolic Serre-Green-Naghdi equations converge (at least formally)
to the original [`SerreGreenNaghdiEquations1D`](@ref). However, the
wave speeds of the hyperbolic system increase with increasing ``\lambda``,
so that explicit time integration methods become more expensive.

Two types of `bathymetry_type` are supported:
- [`bathymetry_flat`](@ref): flat bathymetry (typically ``b = 0`` everywhere)
- [`bathymetry_mild_slope`](@ref): variable bathymetry with mild-slope approximation

For the mild-slope approximation, the Serre-Green-Naghdi equations are
```math
\begin{aligned}
  h_t + (h v)_x &= 0,\\
  h v_t + \frac{1}{2} g (h^2)_x + \frac{1}{2} h (v^2)_x
    + \biggl( \frac{\lambda}{3} H (1 - H / h) \biggr)_x
    + \biggl( g h + \frac{\lambda}{2} (1 - H / h) \biggr) b_x &= 0,\\
  h w_t + h v w_x &= \lambda (1 - H / h),\\
  H_t + H_x u + \frac{3}{2} b_x v &= w.
\end{aligned}
```

References for the hyperbolized Serre-Green-Naghdi system can be found in
- Favrie and Gavrilyuk.
  A rapid numerical method for solving Serre-Green-Naghdi equations
  describing long free surface gravity waves
  [DOI: 10.1088/1361-6544/aa712d](https://doi.org/10.1088/1361-6544/aa712d)
- Busto, Dumbser, Escalante, Favrie, and Gavrilyuk.
  On High Order ADER Discontinuous Galerkin Schemes for First Order Hyperbolic
  Reformulations of Nonlinear Dispersive Systems
  [DOI: 10.1007/s10915-021-01429-8](https://doi.org/10.1007/s10915-021-01429-8)

The semidiscretization implemented here conserves
- the total water mass (integral of ``h``) as a linear invariant
- the total modified energy

for periodic boundary conditions (see Ranocha and Ricchiuto (2024)).
Additionally, it is well-balanced for the lake-at-rest stationary solution, see
- Hendrik Ranocha and Mario Ricchiuto (2024)
  Structure-preserving approximations of the Serre-Green-Naghdi
  equations in standard and hyperbolic form
  [arXiv: 2408.02665](https://arxiv.org/abs/2408.02665)
"""
struct HyperbolicSerreGreenNaghdiEquations1D{Bathymetry <:
                                             Union{BathymetryFlat, BathymetryMildSlope},
                                             RealT <: Real} <:
       AbstractSerreGreenNaghdiEquations{1, 5}
    bathymetry_type::Bathymetry # type of bathymetry
    gravity::RealT # gravitational acceleration
    eta0::RealT # constant still-water surface
    lambda::RealT # hyperbolic relaxation parameter (→ ∞ for Serre-Green-Naghdi)
end

function HyperbolicSerreGreenNaghdiEquations1D(; bathymetry_type = bathymetry_mild_slope,
                                               gravity,
                                               eta0 = 0.0,
                                               lambda)
    HyperbolicSerreGreenNaghdiEquations1D(bathymetry_type, gravity, eta0, lambda)
end

function varnames(::typeof(prim2prim), ::HyperbolicSerreGreenNaghdiEquations1D)
    return ("η", "v", "D", "w", "H")
end
function varnames(::typeof(prim2cons), ::HyperbolicSerreGreenNaghdiEquations1D)
    return ("h", "hv", "b", "hw", "hH")
end

"""
    prim2phys(q, equations::HyperbolicSerreGreenNaghdiEquations1D)

Return the physical variables ``\\eta, v, D`` used also by the
[`SerreGreenNaghdiEquations1D`](@ref) from the main variables
`q` for the [`HyperbolicSerreGreenNaghdiEquations1D`](@ref).
"""
function prim2phys(q, ::HyperbolicSerreGreenNaghdiEquations1D)
    eta, v, D = q
    return SVector(eta, v, D)
end

varnames(::typeof(prim2phys), ::HyperbolicSerreGreenNaghdiEquations1D) = ("η", "v", "D")

is_hyperbolic_appproximation(::HyperbolicSerreGreenNaghdiEquations1D) = Val{true}()

function hyperbolic_approximation_limit(equations::HyperbolicSerreGreenNaghdiEquations1D)
    (; bathymetry_type, gravity, eta0) = equations
    return SerreGreenNaghdiEquations1D(bathymetry_type, gravity, eta0)
end

# compute the auxiliary variables from the main ones
# using the default initialization
function set_approximation_variables!(q, mesh,
                                      equations::HyperbolicSerreGreenNaghdiEquations1D,
                                      solver)
    (; D1) = solver
    eta, v, D, w, H = q.x

    # H ≈ h
    @.. H = eta + D - equations.eta0 # h = (h + b) + (eta0 - b) - eta0

    # w ≈ -h v_x + 3/2 v b_x
    mul!(w, D1, v)
    @.. w = -H * w
    if !(equations.bathymetry_type isa BathymetryFlat)
        b = eta - H
        b_x = D1 * b
        @.. w += 1.5 * v * b_x
    end

    return nothing
end

# TODO: There is a name clash. For the SerreGreenNaghdiEquations1D,
#       the corresponding function is called initial_condition_convergence_test
#       However, we cannot use that name since it's not an analytical solution.
#       How shall we handle this?
"""
    initial_condition_soliton(x, t, equations::HyperbolicSerreGreenNaghdiEquations1D, mesh)

A soliton solution of the [`SerreGreenNaghdiEquations1D`](@ref)
used for convergence tests in a periodic domain. This is physically the
same as [`initial_condition_convergence_test`](@ref) for the
[`SerreGreenNaghdiEquations1D`](@ref). Please note that this is not an
exact solution of the [`HyperbolicSerreGreenNaghdiEquations1D`](@ref)
(only in the limit of the relaxation parameter ``\\lambda \\to \\infty``).

See also [`initial_condition_convergence_test`](@ref).
"""
function initial_condition_soliton(x, t, equations::HyperbolicSerreGreenNaghdiEquations1D,
                                   mesh)
    g = gravity(equations)

    # setup parameters data
    h1 = 1.0
    h2 = 1.2
    c = sqrt(g * h2)

    x_t = mod(x - c * t - xmin(mesh), xmax(mesh) - xmin(mesh)) + xmin(mesh)

    h = h1 + (h2 - h1) * sech(x_t / 2 * sqrt(3 * (h2 - h1) / (h1^2 * h2)))^2
    v = c * (1 - h1 / h)
    # w = -h v_x (b = 0)
    w = -h1 * sqrt(g * h2) * sqrt((-3 * h1 + 3 * h2) / (h1^2 * h2)) * (-h1 + h2) *
        (-h1 - (-h1 + h2) * sech(x_t * sqrt((-3 * h1 + 3 * h2) / (h1^2 * h2)) / 2)^2) *
        tanh(x_t * sqrt((-3 * h1 + 3 * h2) / (h1^2 * h2)) / 2) *
        sech(x_t * sqrt((-3 * h1 + 3 * h2) / (h1^2 * h2)) / 2)^2 /
        (h1 + (-h1 + h2) * sech(x_t * sqrt((-3 * h1 + 3 * h2) / (h1^2 * h2)) / 2)^2)^2
    H = h

    return SVector(h, v, 0, w, H)
end

"""
    initial_condition_manufactured(x, t, equations::HyperbolicSerreGreenNaghdiEquations1D, mesh)

A smooth manufactured solution in combination with
[`source_terms_manufactured`](@ref), see
- Hendrik Ranocha and Mario Ricchiuto (2024)
  Structure-preserving approximations of the Serre-Green-Naghdi
  equations in standard and hyperbolic form
  [arXiv: 2408.02665](https://arxiv.org/abs/2408.02665)
"""
function initial_condition_manufactured(x, t,
                                        equations::HyperbolicSerreGreenNaghdiEquations1D,
                                        mesh)
    eta = 2 + cospi(2 * (x - 2 * t))
    b = -5 - 2 * cospi(2 * x)
    h = eta - b
    v = sinpi(2 * (x - t / 2))
    D = equations.eta0 - b
    # w = -h v_x + 3/2 v b_x
    w = -h * 2 * pi * cospi(2 * (x - t / 2)) + 3 / 2 * v * (4 * pi * sinpi(2 * x))
    H = h
    return SVector(eta, v, D, w, H)
end

"""
    source_terms_manufactured(q, x, t, equations::HyperbolicSerreGreenNaghdiEquations1D)

A smooth manufactured solution in combination with
[`initial_condition_manufactured`](@ref).
"""
function source_terms_manufactured(q, x, t,
                                   equations::HyperbolicSerreGreenNaghdiEquations1D)
    g = gravity(equations)

    a1 = sinpi(4 * t - 2 * x)
    a2 = cospi(4 * t - 2 * x)
    a5 = sinpi(t - 2 * x)
    a6 = cospi(t - 2 * x)
    a7 = sinpi(2 * x)
    a8 = cospi(2 * x)
    a9 = (-a2 - 2 * a8 - 7)
    a10 = (2 * pi * a1 - 4 * pi * a7)

    # Source terms for variable bathymetry
    s1 = -4 * pi * a1 - a5 * a10 - 2 * pi * a6 * a9
    s2 = -2 * pi * a5 * a6 - pi * a6 + 4 * pi * a7 * g + g * a10
    s3 = zero(s1)
    s4 = 8 * pi^2 * a1 * a6 - 2 * pi^2 * a5 * a9 -
         a5 * (-12.0 * pi^2 * a5 * a8 + 4 * pi^2 * a5 * a9 +
          12.0 * pi^2 * a6 * a7 - 2 * pi * a6 * a10) -
         6.0 * pi^2 * a6 * a7
    s5 = -4 * pi * a1 - a5 * a10 - 2 * pi * a6 * a9

    return SVector(s1, s2, s3, s4, s5)
end

"""
    initial_condition_manufactured_reflecting(x, t, equations::HyperbolicSerreGreenNaghdiEquations1D, mesh)

A smooth manufactured solution for reflecting boundary conditions in combination
with [`source_terms_manufactured_reflecting`](@ref).
"""
function initial_condition_manufactured_reflecting(x, t,
                                                   equations::HyperbolicSerreGreenNaghdiEquations1D,
                                                   mesh)
    h = 1 + (2 * t) * (cospi(x) + x + 2)
    v = (-t * x) * sinpi(x)
    D = zero(h) # flat bathymetry (D = 0)
    if equations.bathymetry_type isa BathymetryMildSlope
        D = -(2 * x)
    end
    # w = -h v_x + 3/2 v b_x
    w = (-t * sinpi(x) - pi * t * x * cospi(x)) * (-1 - 2 * t * (2 + x + cospi(x)))
    H = h

    b = equations.eta0 - D
    eta = h + b

    return SVector(eta, v, D, w, H)
end

"""
    source_terms_manufactured_reflecting(q, x, t, equations::HyperbolicSerreGreenNaghdiEquations1D)

A smooth manufactured solution in combination with
[`initial_condition_manufactured_reflecting`](@ref).

calculated as shown in:
https://github.com/NumericalMathematics/DispersiveShallowWater.jl/pull/228#issuecomment-3123350726
"""
function source_terms_manufactured_reflecting(q, x, t,
                                              equations::HyperbolicSerreGreenNaghdiEquations1D)
    g = gravity(equations)

    a1 = 2 + x + cospi(x)
    a2 = pi * t * x * cospi(x)
    t_2 = t^2
    pi_2 = pi^2

    s1 = 2 * a1 + t * (-1 - 2t * a1) * sinpi(x) +
         pi * t * (-1 - 2t * a1) * x * cospi(x) -
         2(1 - pi * sinpi(x)) * t_2 * x * sinpi(x)
    s2 = 2 * g * t - x * sinpi(x) - 2 * g * pi * t * sinpi(x) +
         t_2 * x * (sinpi(x)^2) + pi * t_2 * (x^2) * cospi(x) * sinpi(x)
    s3 = zero(s1)
    s4 = -2(-t * sinpi(x) - a2) * a1 +
         (-sinpi(x) - pi * x * cospi(x)) * (-1 - 2t * a1) +
         ((2pi * t * cospi(x) - pi_2 * t * x * sinpi(x)) *
          (-1 - 2t * a1) +
          2(-t * sinpi(x) - a2) * (1 - pi * sinpi(x)) * t) * t *
         x * sinpi(x)
    s5 = 2 * a1 +
         (t * sinpi(x) + a2) * (-1 - 2t * a1) -
         2(1 - pi * sinpi(x)) * t_2 * x * sinpi(x)
    if equations.bathymetry_type isa BathymetryMildSlope
        s2 += 2 * g
        s4 -= 3x * sinpi(x) +
              ((2pi * t * cospi(x) - pi_2 * t * x * sinpi(x)) *
               (-1 - 2t * a1) +
               2(-t * sinpi(x) - a2) * (1 - pi * sinpi(x)) * t) * t * x *
              sinpi(x) +
              (-3t * sinpi(x) - 3 * a2 +
               (-2pi * t * cospi(x) + pi_2 * t * x * sinpi(x)) *
               (-1 - 2t * a1) -
               2(-t * sinpi(x) - a2) * (1 - pi * sinpi(x)) * t) * t * x *
              sinpi(x)
    end
    return SVector(s1, s2, s3, s4, s5)
end

function create_cache(mesh, equations::HyperbolicSerreGreenNaghdiEquations1D,
                      solver, initial_condition,
                      boundary_conditions::Union{BoundaryConditionPeriodic,
                                                 BoundaryConditionReflecting},
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
    H_over_h = DiffCache(zero(template), N)
    h_x = DiffCache(zero(template), N)
    v_x = DiffCache(zero(template), N)
    hv_x = DiffCache(zero(template), N)
    v2_x = DiffCache(zero(template), N)
    h_hpb_x = DiffCache(zero(template), N)
    H_x = DiffCache(zero(template), N)
    H2_h_x = DiffCache(zero(template), N)
    w_x = DiffCache(zero(template), N)
    hvw_x = DiffCache(zero(template), N)
    tmp = DiffCache(zero(template), N)

    cache = (; h, b, b_x, H_over_h, h_x, v_x, hv_x, v2_x, h_hpb_x, H_x, H2_h_x, w_x, hvw_x,
             tmp)
    return cache
end

# Discretization that conserves
# - the total water mass (integral of ``h``) as a linear invariant
# - the total modified energy
# for periodic boundary conditions, see
# - Hendrik Ranocha and Mario Ricchiuto (2024)
#   Structure-preserving approximations of the Serre-Green-Naghdi
#   equations in standard and hyperbolic form
#   [arXiv: 2408.02665](https://arxiv.org/abs/2408.02665)
# for reflecting boundary conditions, calculation not published yet.
function rhs!(dq, q, t, mesh,
              equations::HyperbolicSerreGreenNaghdiEquations1D,
              initial_condition,
              boundary_conditions::Union{BoundaryConditionPeriodic,
                                         BoundaryConditionReflecting},
              source_terms,
              solver, cache)
    # Unpack physical parameters and SBP operator `D1`
    g = gravity(equations)
    (; lambda, bathymetry_type) = equations
    (; D1) = solver

    # `q` and `dq` are `ArrayPartition`s. They collect the individual
    # arrays for the total water height `eta = h + b`, the velocity `v`,
    # and the additional variables `w` and `H`.
    eta, v, D, w, H = q.x
    dh, dv, dD, dw, dH = dq.x # dh = deta since b is constant in time
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
        H_over_h = get_tmp(cache.H_over_h, eta)
        h_x = get_tmp(cache.h_x, eta)
        v_x = get_tmp(cache.v_x, eta)
        hv_x = get_tmp(cache.hv_x, eta)
        v2_x = get_tmp(cache.v2_x, eta)
        h_hpb_x = get_tmp(cache.h_hpb_x, eta)
        H_x = get_tmp(cache.H_x, eta)
        H2_h_x = get_tmp(cache.H2_h_x, eta)
        w_x = get_tmp(cache.w_x, eta)
        hvw_x = get_tmp(cache.hvw_x, eta)
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

        # hv2_x = D1 * (h * v)
        @.. tmp = h * v
        mul!(hv_x, D1, tmp)

        # v2_x = D1 * (v.^2)
        @.. tmp = v^2
        mul!(v2_x, D1, tmp)

        # h_hpb_x = D1 * (h .* eta)
        @.. tmp = h * eta
        mul!(h_hpb_x, D1, tmp)

        # H_x = D1 * H
        mul!(H_x, D1, H)

        # H2_h_x = D1 * (H^2 / h)
        @.. H_over_h = H / h
        @.. tmp = H * H_over_h
        mul!(H2_h_x, D1, tmp)

        # w_x = D1 * w
        mul!(w_x, D1, w)

        # hvw_x = D1 * (h * v * w)
        @.. tmp = h * v * w
        mul!(hvw_x, D1, tmp)

        # Plain: h_t + (h v)_x = 0
        #
        # Split form for energy conservation:
        # h_t + h_x v + h v_x = 0
        @.. dh = -(h_x * v + h * v_x)
        if boundary_conditions isa BoundaryConditionReflecting
            dh[1] -= h[1] * v[1] / left_boundary_weight(D1)
            dh[end] += h[end] * v[end] / right_boundary_weight(D1)
        end

        # Plain: h v_t + h v v_x + g (h + b) h_x
        #              + ... = 0
        #
        # Split form for energy conservation:
        # h v_t + g (h (h + b))_x - g (h + b) h_x
        #       + 1/2 h (v^2)_x - 1/2 v^2 h_x  + 1/2 v (h v)_x - 1/2 h v v_x
        #       + λ/6 H^2 / h^2 h_x + λ/3 H_x - λ/3 H/h H_x - λ/6 (H^2 / h)_x
        #       + λ/2 b_x - λ/2 H/h b_x = 0
        lambda_6 = lambda / 6
        lambda_3 = lambda / 3
        @.. dv = -(g * h_hpb_x - g * eta * h_x
                   + 0.5 * (h * v2_x - v^2 * h_x)
                   + 0.5 * v * (hv_x - h * v_x)
                   + lambda_6 * (H_over_h^2 * h_x - H2_h_x)
                   + lambda_3 * (1 - H_over_h) * H_x) / h
        if !(bathymetry_type isa BathymetryFlat)
            lambda_2 = lambda / 2
            @.. dv -= lambda_2 * (1 - H_over_h) * b_x / h
        end

        # Plain: h w_t + h v w_x = λ - λ H / h
        #
        # Split form for energy conservation:
        # h w_t + 1/2 (h v w)_x + 1/2 h v w_x
        #       - 1/2 h_x v w - 1/2 h w v_x = λ - λ H / h
        @.. dw = (-0.5 * (hvw_x + h * v * w_x + dh * w) +
                  lambda * (1 - H_over_h)) / h

        # No special split form for energy conservation required:
        # H_t + v H_x + 3/2 v b_x = w
        @.. dH = w - v * H_x
        if !(bathymetry_type isa BathymetryFlat)
            @.. dH -= 1.5 * v * b_x
        end
    end

    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)

    return nothing
end

@inline function prim2cons(q, equations::HyperbolicSerreGreenNaghdiEquations1D)
    h = waterheight(q, equations)
    v = velocity(q, equations)
    b = bathymetry(q, equations)

    hv = h * v
    hw = h * q[4]
    hH = h * q[5]
    return SVector(h, hv, b, hw, hH)
end

@inline function cons2prim(u, equations::HyperbolicSerreGreenNaghdiEquations1D)
    h, hv, b, hw, hH = u

    eta = h + b
    v = hv / h
    D = equations.eta0 - b
    w = hw / h
    H = hH / h
    return SVector(eta, v, D, w, H)
end

# The entropy/energy takes the whole `q` for every point in space
"""
    DispersiveShallowWater.energy_total_modified!(e, q_global, equations::HyperbolicSerreGreenNaghdiEquations1D, cache)

Return the modified total energy `e` of the primitive variables `q_global` for the
[`HyperbolicSerreGreenNaghdiEquations1D`](@ref).
It contains additional terms compared to the usual [`energy_total`](@ref)
modeling non-hydrostatic contributions. The `energy_total_modified`
is a conserved quantity (for periodic boundary conditions).

For a [`bathymetry_mild_slope`](@ref) (and a [`bathymetry_flat`](@ref)),
the total modified energy is given by
```math
\\frac{1}{2} g \\eta^2 + \\frac{1}{2} h v^2 +
\\frac{1}{6} h w^2 + \\frac{\\lambda}{6} h (1 - \\eta / h)^2.
```

`q_global` is a vector of the primitive variables at ALL nodes.

See also [`energy_total_modified`](@ref).
"""
function energy_total_modified!(e, q_global,
                                equations::HyperbolicSerreGreenNaghdiEquations1D,
                                cache)
    # unpack physical parameters and SBP operator `D1`
    g = gravity(equations)
    (; lambda) = equations
    h = get_tmp(cache.h, q_global)
    b = get_tmp(cache.b, q_global)

    # `q_global` is an `ArrayPartition`. It collects the individual arrays for
    # the total water height `eta = h + b` and the velocity `v`.
    eta, v, D, w, H = q_global.x
    @.. b = equations.eta0 - D
    @.. h = eta - b

    # 1/2 g eta^2 + 1/2 h v^2 + 1/6 h^3 w^2 + λ/6 h (1 - H/h)^2
    @.. e = 1 / 2 * g * eta^2 + 1 / 2 * h * v^2 + 1 / 6 * h * w^2 +
            lambda / 6 * h * (1 - H / h)^2

    return e
end
