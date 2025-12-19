"""
    RelaxationCallback(invariant)

Use a relaxation method in time in order to exactly preserve the (nonlinear)
`invariant` for a conservative semidiscretization. A possible choice for
`invariant` is `invariant = entropy`.

Reference
- Hendrik Ranocha, Mohammed Sayyari, Lisandro Dalcin, Matteo Parsani, David I. Ketcheson (2020)
  Relaxation Runge–Kutta Methods: Fully-Discrete Explicit Entropy-Stable Schemes for the
  Compressible Euler and Navier–Stokes Equations
  [DOI: 10.1137/19M1263480](https://doi.org/10.1137/19M1263480)
"""
mutable struct RelaxationCallback{Invariant}
    invariant::Invariant
end

function Base.show(io::IO, cb::DiscreteCallback{<:Any, <:RelaxationCallback})
    @nospecialize cb # reduce precompilation time

    relaxation_callback = cb.affect!
    @unpack invariant = relaxation_callback
    print(io, "RelaxationCallback(invariant=", string(nameof(invariant)), ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain",
                   cb::DiscreteCallback{<:Any, <:RelaxationCallback})
    @nospecialize cb # reduce precompilation time

    if get(io, :compact, false)
        show(io, cb)
    else
        relaxation_callback = cb.affect!

        println(io, "RelaxationCallback")
        print(io, "    invariant: ", string(nameof(relaxation_callback.invariant)))
    end
end

function RelaxationCallback(; invariant)
    relaxation_callback = RelaxationCallback(invariant)

    return DiscreteCallback(relaxation_callback, relaxation_callback, # the first one is the condition, the second the affect!
                            save_positions = (false, false),
                            initialize = initialize!)
end

function initialize!(cb::DiscreteCallback{Condition, Affect!}, u, t,
                     integrator) where {Condition, Affect! <: RelaxationCallback}
    return nothing
end

# this method is called to determine whether the callback should be activated
function (relaxation_callback::RelaxationCallback)(u, t, integrator)
    return true
end

# This method is called as callback during the time integration.
@inline function (relaxation_callback::RelaxationCallback)(integrator)
    semi = integrator.p
    told = integrator.tprev
    qold = integrator.uprev
    tnew = integrator.t
    qnew = integrator.u

    terminate_integration = false
    gamma_lo = one(tnew) / 2
    gamma_hi = 3 * one(tnew) / 2

    @unpack tmp1 = semi.cache # of size N
    @unpack tmp_partitioned = semi.cache # of size nvariables * N and ArrayPartition

    function relaxation_functional(tmp1, q, semi)
        return integrate_quantity!(tmp1, relaxation_callback.invariant, q, semi)
    end

    function convex_combination!(tmp_partitioned, gamma, uold, unew) # for arrays
        @.. tmp_partitioned = uold + gamma * (unew - uold)
        return nothing
    end

    energy_old = relaxation_functional(tmp1, qold, semi)

    # define the root function
    function root(g, tmp1, tmp_partitioned, qold, qnew, energy_old)
        convex_combination!(tmp_partitioned, g, qold, qnew)
        return (relaxation_functional(tmp1, tmp_partitioned, semi) - energy_old)
    end
    # close it over the parameters for less allocations
    root_closure(g) = root(g, tmp1, tmp_partitioned, qold, qnew, energy_old)

    @trixi_timeit timer() "relaxation" begin
        convex_combination!(tmp_partitioned, gamma_lo, qold, qnew)
        val1 = relaxation_functional(tmp1, tmp_partitioned, semi) - energy_old

        convex_combination!(tmp_partitioned, gamma_hi, qold, qnew)
        val2 = relaxation_functional(tmp1, tmp_partitioned, semi) - energy_old

        if (val1 * val2) > 0
            # terminate if there is no root in [gamma_lo, gamma_hi]
            gamma = zero(typeof(tnew))
        else
            gamma = find_zero(root_closure, (gamma_lo, gamma_hi), AlefeldPotraShi())
        end

        if gamma < eps(typeof(gamma))
            terminate_integration = true
        end

        convex_combination!(qnew, gamma, qold, qnew)
        DiffEqBase.set_u!(integrator, qnew)

        if !isapprox(tnew, first(integrator.opts.tstops))
            tgamma = told + gamma * (tnew - told) # convex_combination scalar eq.
            DiffEqBase.set_t!(integrator, tgamma)
        end

        if terminate_integration
            terminate!(integrator)
        end
    end
    return nothing
end
