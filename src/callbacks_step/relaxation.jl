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

    DiscreteCallback(relaxation_callback, relaxation_callback, # the first one is the condition, the second the affect!
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
    tmp222 = similar(qold) # of size nvariables * N
    # @unpack tmp222 = semi.cache # of size nvariables * N and ArrayPartition


    function relaxation_functional(tmp1, q, semi)
        return integrate_quantity!(tmp1, relaxation_callback.invariant, q, semi)
    end

    function convex_combination!(tmp222, gamma, old, new)
        @.. tmp222 = old + gamma * (new - old)
        return nothing
    end

    function convex_combination(gamma, old, new)
        return @.. old + gamma * (new - old)
    end

    function root(g)
        convex_combination!(tmp222, g, qold, qnew)
        return (relaxation_functional(tmp1, tmp222 ,semi) - energy_old)
    end

    energy_old = relaxation_functional(tmp1, qold, semi)

    @trixi_timeit timer() "relaxation" begin
        convex_combination!(tmp222, gamma_lo, qold, qnew)
        val1 = relaxation_functional(tmp1, tmp222, semi) - energy_old

        convex_combination!(tmp222, gamma_hi, qold, qnew)
        val2 = relaxation_functional(tmp1, tmp222, semi) - energy_old

        if (val1 * val2) > 0
            terminate_integration = true
        else
            gamma = find_zero(root, (gamma_lo, gamma_hi), AlefeldPotraShi())
        end

        if gamma < eps(typeof(gamma))
            terminate_integration = true
        end

 
        convex_combination!(tmp222, gamma, qold, qnew)
        DiffEqBase.set_u!(integrator, tmp222)

        if !isapprox(tnew, first(integrator.opts.tstops))

            # convex_combination!(tmp222, gamma, told, tnew)
            
            tgamma = convex_combination(gamma, told, tnew)
            DiffEqBase.set_t!(integrator, tgamma)
        end

        if terminate_integration
            terminate!(integrator)
        end
    end
    return nothing
end



#= This method is called as callback during the time integration.
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
    @unpack tmp222 = semi.cache # of size nvariables * N
    tmp2 = tmp222
    tmp3 = similar(tmp2)



    function relaxation_functional(tmp1, q, semi)
        return integrate_quantity!(tmp1, relaxation_callback.invariant, q, semi)
    end

    function convex_combination!(tmp2, gamma, old, new)
        @. tmp2 = old + gamma * (new - old)
        return nothing
    end

    function convex_combination(gamma, old, new)
        return @.. old + gamma * (new - old)
    end

    function root(g)
        convex_combination!(tmp2, g, qold, qnew)
        return (relaxation_functional(tmp1, tmp2 ,semi) - energy_old)
    end

    function relaxation_functional2(q, semi)
        @unpack tmp1 = semi.cache
        return integrate_quantity!(tmp1, relaxation_callback.invariant, q, semi)
    end

    energy_old = relaxation_functional(tmp1, qold, semi)

    @trixi_timeit timer() "relaxation" begin

        convex_combination!(tmp2, gamma_hi, qold, qnew)
        val2 = relaxation_functional2(tmp2, semi) - energy_old

        convex_combination!(tmp3, gamma_lo, qold, qnew)
        val1 = relaxation_functional2(tmp3, semi) - energy_old

        
        teststs = 
        val1_dif = relaxation_functional2(convex_combination(gamma_lo, qold, qnew), semi) - energy_old
        val2_dif = relaxation_functional2(convex_combination(gamma_hi, qold, qnew), semi) -energy_old

        @show val1
        @show val1_dif
        @show val2
        @show val2_dif


        if (val1 * val2) > 0
            terminate_integration = true
        else
            gamma = find_zero(root, (gamma_lo, gamma_hi), AlefeldPotraShi())
            @show gamma
        end

        if gamma < eps(typeof(gamma))
            terminate_integration = true
        end

 
        convex_combination!(tmp2, gamma, qold, qnew)
        DiffEqBase.set_u!(integrator, tmp2)

        if !isapprox(tnew, first(integrator.opts.tstops))

            # convex_combination!(tmp2, gamma, told, tnew)

            tgamma = convex_combination(gamma, told, tnew)
            DiffEqBase.set_t!(integrator, tgamma)
        end

        if terminate_integration
            terminate!(integrator)
        end
    end
    @show "einmal durch"
    return nothing
end
=#