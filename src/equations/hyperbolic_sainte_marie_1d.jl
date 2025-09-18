struct HyperbolicSainteMarieEquations1D{Bathymetry <: Union{BathymetryFlat, BathymetryMildSlope},
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
                                            alpha)
    c_squared = alpha^2 * gravity * eta0
    HyperbolicSainteMarieEquations1D(bathymetry_type, gravity, eta0, alpha, c_squared)
end

# TODO: What is easier - hv or v as main computational variable?
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
    @.. p = 0

    # w ≈ -(1/2) h v_x + v b_x
    # with h = (h + b) + (eta0 - b) - eta0
    mul!(w, D1, v)
    @.. w = -0.5 * (eta + D - equations.eta0) * w
    if !(equations.bathymetry_type isa BathymetryFlat)
        b = equations.eta0 .- D
        b_x = D1 * b
        @.. w += 1.5 * v * b_x
    end

    return nothing
end

# TODO: other methods

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
            error()
        end

        # Plain: h v_t + h v v_x + g (h + b) h_x
        #              + ... = 0
        #
        # Split form for energy conservation:
        # h v_t + g (h (h + b))_x - g (h + b) h_x
        #       + 1/2 h (v^2)_x - 1/2 v^2 h_x  + 1/2 v (h v)_x - 1/2 h v v_x
        #       + (h p)_x + 2 p b_x = 0
        # TODO: non-hydrostatic terms for variable bathymetry
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
                  - c_squared * h * v_x - 2 * c_squared * w
                  + 2 * c_squared * v * b_x) / h
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
