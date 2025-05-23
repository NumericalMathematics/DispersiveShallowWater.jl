@doc raw"""
    SerreGreenNaghdiEquations1D(; bathymetry_type = bathymetry_variable,
                                gravity, eta0 = 0.0)

Serre-Green-Naghdi system in one spatial dimension.
The equations for flat bathymetry are given by
```math
\begin{aligned}
  h_t + (h v)_x &= 0,\\
  h v_t - \frac{1}{3} (h^3 v_{tx})_x + \frac{1}{2} g (h^2)_x + \frac{1}{2} h (v^2)_x + p_x &= 0,\\
  p &= \frac{1}{3} h^3 v_{x}^2 - \frac{1}{3} h^3 v v_{xx}.
\end{aligned}
```
The unknown quantities of the Serre-Green-Naghdi equations are the
total water height ``\eta = h + b`` and the velocity ``v``.
The gravitational acceleration `gravity` is denoted by ``g`` and the bottom topography
(bathymetry) ``b = \eta_0 - D``. The water height above the bathymetry
is therefore given by ``h = \eta - \eta_0 + D``.
The total water height is therefore given by ``\eta = h + b``.

Three types of `bathymetry_type` are supported:
- [`bathymetry_flat`](@ref): flat bathymetry (typically ``b = 0`` everywhere)
- [`bathymetry_mild_slope`](@ref): variable bathymetry with mild-slope approximation
- [`bathymetry_variable`](@ref): general variable bathymetry

For the mild-slope approximation, the Serre-Green-Naghdi equations are
```math
\begin{aligned}
  h_t + (h v)_x &= 0,\\
  h v_t - \frac{1}{3} (h^3 v_{tx})_x + \frac{1}{2} (h^2 b_x v_t)_x - \frac{1}{2} h^2 b_x v_{tx} + \frac{3}{4} h b_x^2 v_t
    + \frac{1}{2} g (h^2)_x + g h b_x + \frac{1}{2} h (v^2)_x
    + p_x + \frac{3}{2} \frac{p}{h} b_x &= 0,\\
  p &= \frac{1}{3} h^3 v_{x}^2 - \frac{1}{3} h^3 v v_{xx}
    + \frac{1}{2} h^2 v (b_x v)_x.
\end{aligned}
```
For the general case of variable bathymetry without mild-slope
approximation, the Serre-Green-Naghdi equations are
```math
\begin{aligned}
  h_t + (h v)_x &= 0,\\
  h v_t - \frac{1}{3} (h^3 v_{tx})_x + \frac{1}{2} (h^2 b_x v_t)_x - \frac{1}{2} h^2 b_x v_{tx} + h b_x^2 v_t
    + \frac{1}{2} g (h^2)_x + g h b_x + \frac{1}{2} h (v^2)_x
    + p_x + \frac{3}{2} \frac{p}{h} b_x + \psi b_x &= 0,\\
  p &= \frac{1}{3} h^3 v_{x}^2 - \frac{1}{3} h^3 v v_{xx}
    + \frac{1}{2} h^2 v (b_x v)_x,\\
  \psi &= \frac{1}{4} h v (b_x v)_x.
\end{aligned}
```

References for the Serre-Green-Naghdi system can be found in
- Serre (1953)
  Contribution â l'étude des écoulements permanents et variables dans les canaux
  [DOI: 10.1051/lhb/1953034](https://doi.org/10.1051/lhb/1953034)
- Green and Naghdi (1976)
  A derivation of equations for wave propagation in water of variable depth
  [DOI: 10.1017/S0022112076002425](https://doi.org/10.1017/S0022112076002425)

The semidiscretization implemented here conserves
- the total water mass (integral of ``h``) as a linear invariant
- the total momentum (integral of ``h v``) as a nonlinear invariant if the bathymetry is constant
- the total modified energy

for periodic boundary conditions (see Ranocha and Ricchiuto (2024)).
Additionally, it is well-balanced for the lake-at-rest stationary solution, see
- Hendrik Ranocha and Mario Ricchiuto (2024)
  Structure-preserving approximations of the Serre-Green-Naghdi
  equations in standard and hyperbolic form
  [arXiv: 2408.02665](https://arxiv.org/abs/2408.02665)
"""
struct SerreGreenNaghdiEquations1D{Bathymetry <: AbstractBathymetry, RealT <: Real} <:
       AbstractSerreGreenNaghdiEquations{1, 3}
    bathymetry_type::Bathymetry # type of bathymetry
    gravity::RealT # gravitational acceleration
    eta0::RealT # constant still-water surface
end

function SerreGreenNaghdiEquations1D(; bathymetry_type = bathymetry_variable,
                                     gravity, eta0 = 0.0)
    SerreGreenNaghdiEquations1D(bathymetry_type, gravity, eta0)
end

"""
    initial_condition_soliton(x, t, equations::SerreGreenNaghdiEquations1D, mesh)

A classical soliton solution of the Serre-Green-Naghdi equations. This can be used
for convergence tests in a periodic domain, see [`initial_condition_convergence_test`](@ref).
"""
function initial_condition_soliton(x, t, equations::SerreGreenNaghdiEquations1D, mesh)
    g = gravity(equations)

    # setup parameters data
    h1 = 1.0
    h2 = 1.2
    c = sqrt(g * h2)

    x_t = mod(x - c * t - xmin(mesh), xmax(mesh) - xmin(mesh)) + xmin(mesh)

    h = h1 + (h2 - h1) * sech(x_t / 2 * sqrt(3 * (h2 - h1) / (h1^2 * h2)))^2
    v = c * (1 - h1 / h)

    D = zero(h)
    b = -D
    eta = h + b

    return SVector(eta, v, D)
end

"""
    initial_condition_convergence_test(x, t, equations::SerreGreenNaghdiEquations1D, mesh)

A soliton solution used for convergence tests in a periodic domain. Same as
[`initial_condition_soliton`](@ref) for the [`SerreGreenNaghdiEquations1D`](@ref).
"""
function initial_condition_convergence_test(x, t, equations::SerreGreenNaghdiEquations1D,
                                            mesh)
    return initial_condition_soliton(x, t, equations, mesh)
end

"""
    initial_condition_manufactured(x, t, equations::SerreGreenNaghdiEquations1D, mesh)

A smooth manufactured solution in combination with [`source_terms_manufactured`](@ref).
"""
function initial_condition_manufactured(x, t,
                                        equations::SerreGreenNaghdiEquations1D,
                                        mesh)
    h = 5 + cospi(2 * (x - 2 * t))
    v = (1 + sinpi(2 * (x - t / 2)))

    if equations.bathymetry_type isa BathymetryFlat
        D = zero(h)
    else
        D = -(1.5 - cospi(2 * x))
    end

    b = -D
    # watch out: The rhs! calculates using `h`
    # but the initial condition needs to be given as `eta`.
    eta = h + b

    return SVector(eta, v, D)
end
#=
The source terms where calculated using a CAS, here Symbolics.jl
See code: https://github.com/NumericalMathematics/DispersiveShallowWater.jl/pull/180#discussion_r2068562090
For a chosen h and v, in this case
h = 5 + cos(pi * (2 * (x - 2 * t)))
v = (1 + sin(pi * (2 * (x - t / 2))))

...

dh = (h_t + hv_x) # ↓for flat bathymetry↓
dv = (h * v_t - 1//3 * Dx(h^3 * v_tx) + 1//2 * g * Dx(h^2) + 1//2 * h * Dx(v^2) + p_x)

and after this some substitutions to clean everything up.
=#
"""
    source_terms_manufactured(q, x, t, equations::SerreGreenNaghdiEquations1D, mesh)

A smooth manufactured solution in combination with [`initial_condition_manufactured`](@ref).
"""
function source_terms_manufactured(q, x, t,
                                   equations::SerreGreenNaghdiEquations1D{BathymetryFlat})
    g = gravity(equations)

    a3 = cospi(t - 2x)
    a4 = sinpi(t - 2x)
    a5 = sinpi(2t - 4x)
    a6 = sinpi(4t - 2x)
    a7 = 5 + cospi(4t - 2x)
    a8 = a3 * a7^3
    a9 = 1 - a4
    a10 = 8π^3 * a6 * a7^2

    s1 = 2π * (-a6 - a6 * a4 + a3 * a7)

    s2 = (-π * a3 * a7
          + 2g * π * a6 * a7
          + 2π * a3 * a7 * a9
          -
          (1 / 3) * (-12π^3 * a6 * a4 * a7^2
                     +
                     4π^3 * a8)
          + (4 / 3) * π^3 * a5 * a7^3
          + a10 * a3^2 - a10 * a9 * a4
          +
          (8 / 3) * π^3 * a8 * a9)

    return SVector(s1, s2, zero(s1))
end

function source_terms_manufactured(q, x, t,
                                   equations::SerreGreenNaghdiEquations1D{BathymetryVariable})
    g = gravity(equations)

    a1 = cospi(2x)
    a2 = sinpi(2x)
    a3 = cospi(t - 2x)
    a4 = sinpi(t - 2x)
    a5 = sinpi(2t - 4x)
    a6 = sinpi(4t - 2x)
    a7 = cospi(4t - 2x)

    a8 = 5 + a7
    a11 = 1 - a4

    s1 = 2π * (a6 * a11 + a3 * a8 - 2 * a6)

    s2 = π *
         (2 * g * a2 * a7 + 10 * g * a2 + 2 * g * a6 * a7 + 10 * g * a6 + a3 * a7 + 5 * a3 -
          2 * a3 * a4 * a7 - 10 * a3 * a4) +
         π^3 * (8 * a1 * a2 * a4^2 * a7 + 40 * a1 * a2 * a4^2 - 16 * a1 * a2 * a4 * a7 -
          80 * a1 * a2 * a4 + 8 * a1 * a2 * a7 + 40 * a1 * a2 - 12 * a1 * a3 * a4 * a7^2 -
          120 * a1 * a3 * a4 * a7 - 300 * a1 * a3 * a4 + 10 * a1 * a3 * a7^2 +
          100 * a1 * a3 * a7 + 250 * a1 * a3 + 8 * a1 * a4^2 * a6 * a7 +
          40 * a1 * a4^2 * a6 - 16 * a1 * a4 * a6 * a7 - 80 * a1 * a4 * a6 +
          8 * a1 * a6 * a7 + 40 * a1 * a6 - 8 * a2^2 * a3 * a4 * a7 - 40 * a2^2 * a3 * a4 +
          4 * a2^2 * a3 * a7 + 20 * a2^2 * a3 + 8 * a2 * a3^2 * a7^2 + 80 * a2 * a3^2 * a7 +
          200 * a2 * a3^2 - 8 * a2 * a3 * a4 * a6 * a7 - 40 * a2 * a3 * a4 * a6 +
          4 * a2 * a3 * a6 * a7 + 20 * a2 * a3 * a6 - 4 * a2 * a4^2 * a7^2 -
          40 * a2 * a4^2 * a7 - 100 * a2 * a4^2 + 8 * a2 * a4 * a7^2 + 80 * a2 * a4 * a7 +
          200 * a2 * a4 - 4 * a2 * a7^2 - 40 * a2 * a7 - 100 * a2 + 8 * a3^2 * a6 * a7^2 +
          80 * a3^2 * a6 * a7 + 200 * a3^2 * a6 - 8 * a3 * a4 * a7^3 / 3 -
          40 * a3 * a4 * a7^2 - 200 * a3 * a4 * a7 - 1000 * a3 * a4 / 3 +
          4 * a3 * a7^3 / 3 + 20 * a3 * a7^2 + 100 * a3 * a7 + 500 * a3 / 3 +
          8 * a4^2 * a6 * a7^2 + 80 * a4^2 * a6 * a7 + 200 * a4^2 * a6 -
          4 * a4 * a6 * a7^2 - 40 * a4 * a6 * a7 - 100 * a4 * a6 + 4 * a5 * a7^3 / 3 +
          20 * a5 * a7^2 + 100 * a5 * a7 + 500 * a5 / 3)

    return SVector(s1, s2, zero(s1))
end

function source_terms_manufactured(q, x, t,
                                   equations::SerreGreenNaghdiEquations1D{BathymetryMildSlope})
    g = gravity(equations)

    a1 = cospi(2x)
    a2 = sinpi(2x)
    a3 = cospi(t - 2x)
    a4 = sinpi(t - 2x)
    a5 = sinpi(2t - 4x)
    a6 = sinpi(4t - 2x)
    a7 = cospi(4t - 2x)
    a8 = 5 + a7
    a9 = a6 * a7
    a10 = a6 * a7^2
    a11 = 1 - a4
    a12 = a2 * a4
    a13 = a2 * a4^2

    s1 = 2π * (a6 * a11 + a3 * a8 - 2 * a6)

    s2 = π * (2 * a2 * a7 * g + 10 * a2 * g - 2 * a3 * a4 * a7 - 10 * a3 * a4 + a3 * a7 +
          5 * a3 + 2 * a9 * g + 10 * a6 * g) +
         π^3 * (6 * a1 * a13 * a7 + 30 * a1 * a13 - 12 * a1 * a12 * a7 -
          60 * a1 * a12 + 6 * a1 * a2 * a7 + 30 * a1 * a2 - 12 * a1 * a3 * a4 * a7^2 -
          120 * a1 * a3 * a4 * a7 - 300 * a1 * a3 * a4 + 10 * a1 * a3 * a7^2 +
          100 * a1 * a3 * a7 + 250 * a1 * a3 + 8 * a1 * a4^2 * a9 +
          40 * a1 * a4^2 * a6 - 16 * a1 * a4 * a9 - 80 * a1 * a4 * a6 +
          8 * a1 * a9 + 40 * a1 * a6 - 6 * a2^2 * a3 * a4 * a7 - 30 * a2^2 * a3 * a4 +
          3 * a2^2 * a3 * a7 + 15 * a2^2 * a3 + 8 * a2 * a3^2 * a7^2 + 80 * a2 * a3^2 * a7 +
          200 * a2 * a3^2 - 8 * a2 * a3 * a4 * a9 - 40 * a2 * a3 * a4 * a6 +
          4 * a2 * a3 * a9 + 20 * a2 * a3 * a6 - 4 * a13 * a7^2 -
          40 * a13 * a7 - 100 * a13 + 8 * a12 * a7^2 + 80 * a12 * a7 +
          200 * a12 - 4 * a2 * a7^2 - 40 * a2 * a7 - 100 * a2 + 8 * a3^2 * a10 +
          80 * a3^2 * a9 + 200 * a3^2 * a6 - 8 * a3 * a4 * a7^3 / 3 -
          40 * a3 * a4 * a7^2 - 200 * a3 * a4 * a7 - 1000 * a3 * a4 / 3 +
          4 * a3 * a7^3 / 3 + 20 * a3 * a7^2 + 100 * a3 * a7 + 500 * a3 / 3 +
          8 * a4^2 * a10 + 80 * a4^2 * a6 * a7 + 200 * a4^2 * a6 -
          4 * a4 * a10 - 40 * a4 * a6 * a7 - 100 * a4 * a6 + 4 * a5 * a7^3 / 3 +
          20 * a5 * a7^2 + 100 * a5 * a7 + 500 * a5 / 3)

    return SVector(s1, s2, zero(s1))
end

# flat bathymetry
function create_cache(mesh,
                      equations::SerreGreenNaghdiEquations1D{BathymetryFlat},
                      solver,
                      initial_condition,
                      ::BoundaryConditionPeriodic,
                      RealT, uEltype)
    D1 = solver.D1

    # create temporary storage
    h = ones(RealT, nnodes(mesh))
    b = zero(h)
    h_x = zero(h)
    v_x = zero(h)
    h2_x = zero(h)
    hv_x = zero(h)
    v2_x = zero(h)
    h2_v_vx_x = zero(h)
    h_vx_x = zero(h)
    p_x = zero(h)
    tmp = zero(h)
    M_h = zero(h)
    M_h3_3 = zero(h)

    if D1 isa PeriodicUpwindOperators
        v_x_upwind = zero(h)

        D1mat_minus = sparse(D1.minus)

        system_matrix = let cache = (; M_h, M_h3_3)
            assemble_system_matrix!(cache, h,
                                    D1, D1mat_minus, equations)
        end
        factorization = cholesky(system_matrix)

        cache = (; h, b, h_x, v_x, v_x_upwind, h2_x, hv_x, v2_x,
                 h2_v_vx_x, h_vx_x, p_x, tmp,
                 M_h, M_h3_3,
                 D1, D1mat_minus, factorization)
    else
        if D1 isa FourierDerivativeOperator
            D1mat = Matrix(D1)

            cache = (; h, b, h_x, v_x, h2_x, hv_x, v2_x,
                     h2_v_vx_x, h_vx_x, p_x, tmp,
                     M_h, M_h3_3,
                     D1, D1mat)
        else
            D1mat = sparse(D1)

            system_matrix = let cache = (; M_h, M_h3_3)
                assemble_system_matrix!(cache, h,
                                        D1, D1mat, equations)
            end
            factorization = cholesky(system_matrix)

            cache = (; h, b, h_x, v_x, h2_x, hv_x, v2_x,
                     h2_v_vx_x, h_vx_x, p_x, tmp,
                     M_h, M_h3_3,
                     D1, D1mat, factorization)
        end
    end

    return cache
end

# variable bathymetry
function create_cache(mesh,
                      equations::SerreGreenNaghdiEquations1D,
                      solver,
                      initial_condition,
                      ::BoundaryConditionPeriodic,
                      RealT, uEltype)
    D1 = solver.D1

    # create temporary storage
    h = ones(RealT, nnodes(mesh))
    h_x = zero(h)
    b = zero(h)
    b_x = zero(h)
    v_x = zero(h)
    h_hpb_x = zero(h)
    hv_x = zero(h)
    v2_x = zero(h)
    h2_v_vx_x = zero(h)
    h_vx_x = zero(h)
    p_h = zero(h)
    p_x = zero(h)
    tmp = zero(h)
    M_h_p_h_bx2 = zero(h)
    M_h3_3 = zero(h)
    M_h2_bx = zero(h)

    # b_x appears in the system matrix, so it must not be zero to get the
    # correct pattern of the factorization. Just setting it to unity does
    # also not work since some terms cancel. Thus, b_x needs to be
    # "variable enough" to get the correct pattern of non-zero entries in
    # the system matrix and the factorization.
    # TODO: This is a hack and should be improved. It would be ideal if we
    #       had access to the initial condition and the exact value of the
    #       bathymetry here.
    let x = grid(D1)
        @.. b_x = x^3
    end

    if D1 isa PeriodicUpwindOperators
        v_x_upwind = zero(h)
        p_0 = zero(h)

        D1mat_minus = sparse(D1.minus)

        system_matrix = let cache = (; M_h_p_h_bx2, M_h3_3, M_h2_bx)
            assemble_system_matrix!(cache, h, b_x,
                                    D1, D1mat_minus, equations)
        end
        factorization = cholesky(system_matrix)

        cache = (; h, h_x, v_x, v_x_upwind, h_hpb_x, b, b_x, hv_x, v2_x,
                 h2_v_vx_x, h_vx_x, p_h, p_0, p_x, tmp,
                 M_h_p_h_bx2, M_h3_3, M_h2_bx,
                 D1, D1mat_minus, factorization)
    else
        if D1 isa FourierDerivativeOperator
            D1mat = Matrix(D1)

            cache = (; h, h_x, v_x, h_hpb_x, b, b_x, hv_x, v2_x,
                     h2_v_vx_x, h_vx_x, p_h, p_x, tmp,
                     M_h_p_h_bx2, M_h3_3, M_h2_bx,
                     D1, D1mat)
        else
            D1mat = sparse(D1)

            system_matrix = let cache = (; M_h_p_h_bx2, M_h3_3, M_h2_bx)
                assemble_system_matrix!(cache, h, b_x,
                                        D1, D1mat, equations)
            end
            factorization = cholesky(system_matrix)

            cache = (; h, h_x, v_x, h_hpb_x, b, b_x, hv_x, v2_x,
                     h2_v_vx_x, h_vx_x, p_h, p_x, tmp,
                     M_h_p_h_bx2, M_h3_3, M_h2_bx,
                     D1, D1mat, factorization)
        end
    end

    if equations.bathymetry_type isa BathymetryVariable
        psi = zero(h)
        cache = (; cache..., psi)
    end

    return cache
end

function assemble_system_matrix!(cache, h, D1, D1mat,
                                 ::SerreGreenNaghdiEquations1D{BathymetryFlat})
    (; M_h, M_h3_3) = cache

    @.. M_h = h
    scale_by_mass_matrix!(M_h, D1)
    @.. M_h3_3 = (1 / 3) * h^3
    scale_by_mass_matrix!(M_h3_3, D1)

    # Floating point errors accumulate a bit and the system matrix
    # is not necessarily perfectly symmetric but only up to
    # round-off errors. We wrap it here to avoid issues with the
    # factorization.
    return Symmetric(Diagonal(M_h) + D1mat' * (Diagonal(M_h3_3) * D1mat))
end

# variable bathymetry
function assemble_system_matrix!(cache, h, b_x, D1, D1mat,
                                 equations::SerreGreenNaghdiEquations1D)
    (; M_h_p_h_bx2, M_h3_3, M_h2_bx) = cache

    if equations.bathymetry_type isa BathymetryMildSlope
        factor = 0.75
    elseif equations.bathymetry_type isa BathymetryVariable
        factor = 1.0
    end
    @.. M_h_p_h_bx2 = h + factor * h * b_x^2
    scale_by_mass_matrix!(M_h_p_h_bx2, D1)
    inv3 = 1 / 3
    @.. M_h3_3 = inv3 * h^3
    scale_by_mass_matrix!(M_h3_3, D1)
    @.. M_h2_bx = 0.5 * h^2 * b_x
    scale_by_mass_matrix!(M_h2_bx, D1)
    # Floating point errors accumulate a bit and the system matrix
    # is not necessarily perfectly symmetric but only up to
    # round-off errors. We wrap it here to avoid issues with the
    # factorization.
    return Symmetric(Diagonal(M_h_p_h_bx2)
                     +
                     D1mat' * (Diagonal(M_h3_3) * D1mat
                               -
                               Diagonal(M_h2_bx))
                     -
                     Diagonal(M_h2_bx) * D1mat)
end

# Discretization that conserves
# - the total water mass (integral of ``h``) as a linear invariant
# - the total momentum (integral of ``h v``) as a nonlinear invariant for flat bathymetry
# - the total modified energy
# for periodic boundary conditions, see
# - Hendrik Ranocha and Mario Ricchiuto (2024)
#   Structure-preserving approximations of the Serre-Green-Naghdi
#   equations in standard and hyperbolic form
#   [arXiv: 2408.02665](https://arxiv.org/abs/2408.02665)
function rhs!(dq, q, t, mesh,
              equations::SerreGreenNaghdiEquations1D,
              initial_condition,
              boundary_conditions::BoundaryConditionPeriodic,
              source_terms,
              solver, cache)
    if cache.D1 isa PeriodicUpwindOperators
        rhs_sgn_upwind!(dq, q, t, equations, source_terms, solver, cache,
                        equations.bathymetry_type,
                        boundary_conditions)
    else
        rhs_sgn_central!(dq, q, t, equations, source_terms, solver, cache,
                         equations.bathymetry_type,
                         boundary_conditions)
    end

    return nothing
end

function rhs_sgn_central!(dq, q, t, equations, source_terms, solver, cache,
                          ::BathymetryFlat, boundary_conditions::BoundaryConditionPeriodic)
    # Unpack physical parameters and SBP operator `D1` as well as the
    # SBP operator in sparse matrix form `D1mat`
    g = gravity(equations)
    (; D1, D1mat) = cache

    # `q` and `dq` are `ArrayPartition`s. They collect the individual
    # arrays for the water height `h` and the velocity `v`.
    eta, v, D = q.x
    dh, dv, dD = dq.x # dh = deta since b is constant in time
    fill!(dD, zero(eltype(dD)))

    @trixi_timeit timer() "hyperbolic terms" begin
        # Compute all derivatives required below
        (; h, b, h_x, v_x, h2_x, hv_x, v2_x, h2_v_vx_x,
        h_vx_x, p_x, tmp) = cache

        @.. b = equations.eta0 - D
        @.. h = eta - b

        mul!(h_x, D1, h)
        mul!(v_x, D1, v)
        @.. tmp = h^2
        mul!(h2_x, D1, tmp)
        @.. tmp = h * v
        mul!(hv_x, D1, tmp)
        @.. tmp = v^2
        mul!(v2_x, D1, tmp)

        @.. tmp = h^2 * v * v_x
        mul!(h2_v_vx_x, D1, tmp)
        @.. tmp = h * v_x
        mul!(h_vx_x, D1, tmp)
        inv6 = 1 / 6
        @.. tmp = (0.5 * h^2 * (h * v_x + h_x * v) * v_x
                   -
                   inv6 * h * h2_v_vx_x
                   -
                   inv6 * h^2 * v * h_vx_x)
        mul!(p_x, D1, tmp)

        # Plain: h_t + (h v)_x = 0
        #
        # Split form for energy conservation:
        # h_t + h_x v + h v_x = 0
        @.. dh = -(h_x * v + h * v_x)

        # Plain: h v_t + ... = 0
        #
        # Split form for energy conservation:
        @.. dv = -(g * h2_x - g * h * h_x
                   +
                   0.5 * h * v2_x
                   -
                   0.5 * v^2 * h_x
                   +
                   0.5 * hv_x * v
                   -
                   0.5 * h * v * v_x
                   +
                   p_x)
    end

    # add source term
    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)

    tmp .= dv

    # The code below is equivalent to
    #   dv .= (Diagonal(h) - D1mat * Diagonal(1/3 .* h.^3) * D1mat) \ tmp
    # but faster since the symbolic factorization is reused.
    @trixi_timeit timer() "assembling elliptic operator" begin
        system_matrix = assemble_system_matrix!(cache, h,
                                                D1, D1mat,
                                                equations)
    end

    @trixi_timeit timer() "solving elliptic system" begin
        solve_system_matrix!(dv, system_matrix,
                             tmp, equations, D1, cache, boundary_conditions)
    end

    return nothing
end

function rhs_sgn_upwind!(dq, q, t, equations, source_terms, solver, cache, ::BathymetryFlat,
                         boundary_conditions::BoundaryConditionPeriodic)
    # Unpack physical parameters and SBP operator `D1` as well as the
    # SBP upwind operator in sparse matrix form `D1mat_minus`
    g = gravity(equations)
    (; D1mat_minus) = cache
    D1_upwind = cache.D1
    D1 = D1_upwind.central

    # `q` and `dq` are `ArrayPartition`s. They collect the individual
    # arrays for the water height `h` and the velocity `v`.
    eta, v, D = q.x
    dh, dv, dD = dq.x # dh = deta since b is constant in time
    fill!(dD, zero(eltype(dD)))

    @trixi_timeit timer() "hyperbolic terms" begin
        # Compute all derivatives required below
        (; h, b, h_x, v_x, v_x_upwind, h2_x, hv_x, v2_x,
        h2_v_vx_x, h_vx_x, p_x, tmp) = cache

        @.. b = equations.eta0 - D
        @.. h = eta - b

        mul!(h_x, D1, h)
        mul!(v_x, D1, v)
        mul!(v_x_upwind, D1_upwind.minus, v)
        @.. tmp = h^2
        mul!(h2_x, D1, tmp)
        @.. tmp = h * v
        mul!(hv_x, D1, tmp)
        @.. tmp = v^2
        mul!(v2_x, D1, tmp)

        @.. tmp = h^2 * v * v_x
        mul!(h2_v_vx_x, D1, tmp)
        @.. tmp = h * v_x
        mul!(h_vx_x, D1, tmp)
        # p_+
        @.. tmp = 0.5 * h^2 * (h * v_x + h_x * v) * v_x_upwind
        mul!(p_x, D1_upwind.plus, tmp)
        # p_0
        minv6 = -1 / 6
        @.. tmp = minv6 * (h * h2_v_vx_x
                           +
                           h^2 * v * h_vx_x)
        mul!(p_x, D1, tmp, 1.0, 1.0)

        # Plain: h_t + (h v)_x = 0
        #
        # Split form for energy conservation:
        # h_t + h_x v + h v_x = 0
        @.. dh = -(h_x * v + h * v_x)

        # Plain: h v_t + ... = 0
        #
        # Split form for energy conservation:
        @.. dv = -(g * h2_x - g * h * h_x
                   +
                   0.5 * h * v2_x
                   -
                   0.5 * v^2 * h_x
                   +
                   0.5 * hv_x * v
                   -
                   0.5 * h * v * v_x
                   +
                   p_x)
    end

    # add source term
    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)

    tmp .= dv
    # The code below is equivalent to
    #   dv .= (Diagonal(h) - D1mat_plus * Diagonal(1/3 .* h.^3) * D1mat_minus) \ tmp
    # but faster since the symbolic factorization is reused.
    @trixi_timeit timer() "assembling elliptic operator" begin
        system_matrix = assemble_system_matrix!(cache, h,
                                                D1, D1mat_minus,
                                                equations)
    end

    @trixi_timeit timer() "solving elliptic system" begin
        solve_system_matrix!(dv, system_matrix,
                             tmp, equations, D1, cache, boundary_conditions)
    end

    return nothing
end

function rhs_sgn_central!(dq, q, t, equations, source_terms, solver, cache,
                          ::Union{BathymetryMildSlope, BathymetryVariable},
                          boundary_conditions::BoundaryConditionPeriodic)
    # Unpack physical parameters and SBP operator `D1` as well as the
    # SBP operator in sparse matrix form `D1mat`
    g = gravity(equations)
    (; D1, D1mat) = cache

    # `q` and `dq` are `ArrayPartition`s. They collect the individual
    # arrays for the water height `h` and the velocity `v`.
    eta, v, D = q.x
    dh, dv, dD = dq.x # dh = deta since b is constant in time
    fill!(dD, zero(eltype(dD)))

    @trixi_timeit timer() "hyperbolic terms" begin
        # Compute all derivatives required below
        (; h, h_x, v_x, h_hpb_x, b, b_x, hv_x, v2_x,
        h2_v_vx_x, h_vx_x, p_h, p_x, tmp) = cache
        if equations.bathymetry_type isa BathymetryVariable
            (; psi) = cache
        end

        @.. b = equations.eta0 - D
        @.. h = eta - b
        mul!(b_x, D1, b)

        mul!(h_x, D1, h)
        mul!(v_x, D1, v)
        @.. tmp = h * eta
        mul!(h_hpb_x, D1, tmp)
        @.. tmp = h * v
        mul!(hv_x, D1, tmp)
        @.. tmp = v^2
        mul!(v2_x, D1, tmp)

        @.. tmp = h^2 * v * v_x
        mul!(h2_v_vx_x, D1, tmp)
        @.. tmp = h * v_x
        mul!(h_vx_x, D1, tmp)
        inv6 = 1 / 6
        @.. p_h = (0.5 * h * (h * v_x + h_x * v) * v_x
                   -
                   inv6 * h2_v_vx_x
                   -
                   inv6 * h * v * h_vx_x)
        @.. tmp = h * b_x * v^2
        mul!(p_x, D1, tmp)
        @.. p_h += 0.25 * p_x
        if equations.bathymetry_type isa BathymetryVariable
            @.. psi = 0.125 * p_x
        end
        @.. tmp = b_x * v
        mul!(p_x, D1, tmp)
        @.. p_h += 0.25 * h * v * p_x
        if equations.bathymetry_type isa BathymetryVariable
            @.. psi += 0.125 * h * v * p_x
        end
        @.. p_h = p_h - 0.25 * (h_x * v + h * v_x) * b_x * v
        if equations.bathymetry_type isa BathymetryVariable
            @.. psi -= 0.125 * (h_x * v + h * v_x) * b_x * v
        end
        @.. tmp = p_h * h
        mul!(p_x, D1, tmp)

        # Plain: h_t + (h v)_x = 0
        #
        # Split form for energy conservation:
        # h_t + h_x v + h v_x = 0
        @.. dh = -(h_x * v + h * v_x)

        # Plain: h v_t + ... = 0
        #
        # Split form for energy conservation:
        @.. dv = -(g * h_hpb_x - g * eta * h_x
                   +
                   0.5 * h * v2_x
                   -
                   0.5 * v^2 * h_x
                   +
                   0.5 * hv_x * v
                   -
                   0.5 * h * v * v_x
                   + p_x
                   + 1.5 * p_h * b_x)
        if equations.bathymetry_type isa BathymetryVariable
            @.. dv = dv - psi * b_x
        end
    end

    # add source term
    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)

    tmp .= dv

    # The code below is equivalent to
    #   dv .= (Diagonal(h .+ factor .* h .* b_x.^2) - D1mat * (Diagonal(1/3 .* h.^3) * D1mat - Diagonal(0.5 .* h.^2 .* b_x) * D1mat) \ tmp
    # but faster since the symbolic factorization is reused.
    @trixi_timeit timer() "assembling elliptic operator" begin
        system_matrix = assemble_system_matrix!(cache, h, b_x,
                                                D1, D1mat,
                                                equations)
    end

    @trixi_timeit timer() "solving elliptic system" begin
        solve_system_matrix!(dv, system_matrix,
                             tmp, equations, D1, cache, boundary_conditions)
    end

    return nothing
end

function rhs_sgn_upwind!(dq, q, t, equations, source_terms, solver, cache,
                         ::Union{BathymetryMildSlope, BathymetryVariable},
                         boundary_conditions::BoundaryConditionPeriodic)
    # Unpack physical parameters and SBP operator `D1` as well as the
    # SBP operator in sparse matrix form `D1mat`
    g = gravity(equations)
    (; D1mat_minus) = cache
    D1_upwind = cache.D1
    D1 = D1_upwind.central

    # `q` and `dq` are `ArrayPartition`s. They collect the individual
    # arrays for the water height `h` and the velocity `v`.
    eta, v, D = q.x
    dh, dv, dD = dq.x # dh = deta since b is constant in time
    fill!(dD, zero(eltype(dD)))

    @trixi_timeit timer() "hyperbolic terms" begin
        # Compute all derivatives required below
        (; h, h_x, v_x, v_x_upwind, h_hpb_x, b, b_x, hv_x, v2_x,
        h2_v_vx_x, h_vx_x, p_h, p_0, p_x, tmp) = cache
        if equations.bathymetry_type isa BathymetryVariable
            (; psi) = cache
        end

        @.. b = equations.eta0 - D
        @.. h = eta - b
        mul!(b_x, D1, b)

        mul!(h_x, D1, h)
        mul!(v_x, D1, v)
        mul!(v_x_upwind, D1_upwind.minus, v)
        @.. tmp = h * eta
        mul!(h_hpb_x, D1, tmp)
        @.. tmp = h * v
        mul!(hv_x, D1, tmp)
        @.. tmp = v^2
        mul!(v2_x, D1, tmp)

        @.. tmp = h^2 * v * v_x
        mul!(h2_v_vx_x, D1, tmp)
        @.. tmp = h * v_x
        mul!(h_vx_x, D1, tmp)
        # p_0
        minv6 = -1 / 6
        @.. p_h = minv6 * (h2_v_vx_x
                           +
                           h * v * h_vx_x)
        @.. tmp = h * b_x * v^2
        mul!(p_x, D1, tmp)
        @.. p_h += 0.25 * p_x
        if equations.bathymetry_type isa BathymetryVariable
            @.. psi = 0.125 * p_x
        end
        @.. tmp = b_x * v
        mul!(p_x, D1, tmp)
        @.. p_h += 0.25 * h * v * p_x
        if equations.bathymetry_type isa BathymetryVariable
            @.. psi += 0.125 * h * v * p_x
        end
        @.. p_0 = p_h * h
        mul!(p_x, D1, p_0)
        # p_+
        @.. tmp = (0.5 * h * (h * v_x + h_x * v) * v_x_upwind
                   -
                   0.25 * (h_x * v + h * v_x) * b_x * v)
        if equations.bathymetry_type isa BathymetryVariable
            @.. psi -= 0.125 * (h_x * v + h * v_x) * b_x * v
        end
        @.. p_h = p_h + tmp
        @.. tmp = tmp * h
        mul!(p_x, D1_upwind.plus, tmp, 1.0, 1.0)

        # Plain: h_t + (h v)_x = 0
        #
        # Split form for energy conservation:
        # h_t + h_x v + h v_x = 0
        @.. dh = -(h_x * v + h * v_x)

        # Plain: h v_t + ... = 0
        #
        # Split form for energy conservation:
        @.. dv = -(g * h_hpb_x - g * eta * h_x
                   +
                   0.5 * h * v2_x
                   -
                   0.5 * v^2 * h_x
                   +
                   0.5 * hv_x * v
                   -
                   0.5 * h * v * v_x
                   + p_x
                   + 1.5 * p_h * b_x)
        if equations.bathymetry_type isa BathymetryVariable
            @.. dv = dv - psi * b_x
        end
    end

    # add source term
    @trixi_timeit timer() "source terms" calc_sources!(dq, q, t, source_terms, equations,
                                                       solver)

    tmp .= dv

    # The code below is equivalent to
    #   dv .= (Diagonal(h .+ factor .* h .* b_x.^2) - D1mat * (Diagonal(1/3 .* h.^3) * D1mat - Diagonal(0.5 .* h.^2 .* b_x)) - Diagonal(0.5 .* h.^2 .* b_x) * D1mat) \ tmp
    # but faster since the symbolic factorization is reused.
    @trixi_timeit timer() "assembling elliptic operator" begin
        system_matrix = assemble_system_matrix!(cache, h, b_x,
                                                D1, D1mat_minus,
                                                equations)
    end

    @trixi_timeit timer() "solving elliptic system" begin
        solve_system_matrix!(dv, system_matrix,
                             tmp, equations, D1, cache, boundary_conditions)
    end

    return nothing
end

# The modified entropy/energy takes the whole `q` for every point in space
"""
    DispersiveShallowWater.energy_total_modified!(e, q_global, equations::SerreGreenNaghdiEquations1D, cache)

Return the modified total energy `e` of the primitive variables `q_global` for the
[`SerreGreenNaghdiEquations1D`](@ref).
It contains an additional term containing a
derivative compared to the usual [`energy_total`](@ref) modeling
non-hydrostatic contributions. The `energy_total_modified`
is a conserved quantity (for periodic boundary conditions).

For a [`bathymetry_flat`](@ref) the total modified energy is given by
```math
\\frac{1}{2} g \\eta^2 + \\frac{1}{2} h v^2 + \\frac{1}{6} h^3 v_x^2.
```
For a [`bathymetry_mild_slope`](@ref) the total modified energy is given by
```math
\\frac{1}{2} g \\eta^2 + \\frac{1}{2} h v^2 + \\frac{1}{6} h (-h v_x + 1.5 v b_x)^2.
```
For a [`bathymetry_variable`](@ref) the total modified energy has the additional term
```math
+ \\frac{1}{8} h (v b_x)^2.
```

`q_global` is a vector of the primitive variables at ALL nodes.
`cache` needs to hold the SBP operators used by the `solver`.

See also [`energy_total_modified`](@ref).
"""
function energy_total_modified!(e, q_global,
                                equations::SerreGreenNaghdiEquations1D,
                                cache)
    # unpack physical parameters and SBP operator `D1`
    g = gravity(equations)
    (; D1, h, b, v_x) = cache

    # `q_global` is an `ArrayPartition`. It collects the individual arrays for
    # the total water height `eta = h + b` and the velocity `v`.
    eta, v, D = q_global.x
    @.. b = equations.eta0 - D
    @.. h = eta - b

    # 1/2 g eta^2 + 1/2 h v^2 + 1/6 h^3 w^2
    # and + 1/8 h (v b_x)^2 for full bathymetry without mild-slope approximation
    if D1 isa PeriodicUpwindOperators
        mul!(v_x, D1.minus, v)
    else
        mul!(v_x, D1, v)
    end

    if equations.bathymetry_type isa BathymetryFlat
        b_x = cache.tmp
        fill!(b_x, zero(eltype(b_x)))
    else
        (; b, b_x) = cache
        @.. b = equations.eta0 - q_global.x[3]
        if D1 isa PeriodicUpwindOperators
            mul!(b_x, D1.central, b)
        else
            mul!(b_x, D1, b)
        end
    end

    @.. e = 1 / 2 * g * eta^2 + 1 / 2 * h * v^2 + 1 / 6 * h * (-h * v_x + 1.5 * v * b_x)^2
    if equations.bathymetry_type isa BathymetryVariable
        @.. e += 1 / 8 * h * (v * b_x)^2
    end

    return e
end
