"""
    AbstractSolver

An abstract supertype of specific solvers.
"""
abstract type AbstractSolver end

"""
    Solver

A `struct` that holds the summation-by-parts (SBP) operators that are used for the spatial discretization.
"""
struct Solver{RealT <: Real,
              FirstDerivative <: AbstractDerivativeOperator{RealT},
              SecondDerivative <:
              Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}, Nothing},
              ThirdDerivative <:
              Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT}, Nothing}} <:
       AbstractSolver
    D1::FirstDerivative
    D2::SecondDerivative
    D3::ThirdDerivative

    function Solver{RealT,
                    FirstDerivative,
                    SecondDerivative,
                    ThirdDerivative}(D1::FirstDerivative,
                                     D2::SecondDerivative,
                                     D3::ThirdDerivative) where {
                                                                 RealT,
                                                                 FirstDerivative,
                                                                 SecondDerivative,
                                                                 ThirdDerivative
                                                                 }
        @assert derivative_order(D1) == 1

        if D2 isa AbstractDerivativeOperator &&
           !(D2 isa SummationByPartsOperators.FourierPolynomialDerivativeOperator) &&
           !(D2 isa SummationByPartsOperators.PeriodicRationalDerivativeOperator)
            @assert derivative_order(D2) == 2
        end

        if D3 isa AbstractDerivativeOperator &&
           !(D3 isa SummationByPartsOperators.FourierPolynomialDerivativeOperator) &&
           !(D3 isa SummationByPartsOperators.PeriodicRationalDerivativeOperator)
            @assert derivative_order(D3) == 3
        end

        new(D1, D2, D3)
    end
end

"""
    Solver(mesh, accuracy_order)

Create a solver, where the three summation-by-parts (SBP) first-, second-, and third-derivative operators
are of accuracy order `accuracy_order` and associated to the `mesh`.

!!! warning "Periodic operators only"
    This constructor creates periodic derivative operators that are only compatible with periodic
    boundary conditions. For non-periodic boundary conditions, use the `Solver(D1, D2, D3)`
    constructor with appropriate non-periodic operators (or `nothing`).
"""
function Solver(mesh, accuracy_order)
    if isodd(accuracy_order)
        @warn "DispersiveShallowWater.jl is expecting a central difference operator. This is not given for an odd accuracy (got $accuracy_order) order.\n This can lead to a significant reduction in the order of convergence of the solution."
    end
    D1 = periodic_derivative_operator(1, accuracy_order, mesh.xmin, mesh.xmax, mesh.N)
    D2 = periodic_derivative_operator(2, accuracy_order, mesh.xmin, mesh.xmax, mesh.N)
    D3 = periodic_derivative_operator(3, accuracy_order + 2, mesh.xmin, mesh.xmax, mesh.N)
    @assert real(D1) == real(D2)
    Solver{real(D1), typeof(D1), typeof(D2), typeof(D3)}(D1, D2, D3)
end

# Also allow to pass custom SBP operators (for convenience without explicitly specifying the type)
"""
    Solver(D1, D2 = nothing, D3 = nothing)

Create a solver, where `D1` is an `AbstractDerivativeOperator`
from [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl)
of first `derivative_order`, `D2` is an `AbstractDerivativeOperator`
of second `derivative_order` or an `AbstractMatrix`, and `D3` is an `AbstractDerivativeOperator`
of third `derivative_order` or an `AbstractMatrix`. `D2` and `D3` can also be `nothing`
if that derivative is not used by the discretization.
All given summation-by-parts operators should be associated with the same grid.
"""
function Solver(D1::AbstractDerivativeOperator{RealT},
                D2::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT},
                          Nothing} = nothing,
                D3::Union{AbstractDerivativeOperator{RealT}, AbstractMatrix{RealT},
                          Nothing} = nothing) where {RealT}
    Solver{RealT, typeof(D1), typeof(D2), typeof(D3)}(D1, D2, D3)
end

function Base.show(io::IO, solver::Solver{RealT}) where {RealT}
    print(io, "Solver{", RealT, "}")
end

function Base.show(io::IO, ::MIME"text/plain", solver::Solver{RealT}) where {RealT}
    if get(io, :compact, false)
        show(io, solver)
    else
        println(io, "Solver{", RealT, "}")
        println(io, "    D1: ", solver.D1)
        println(io, "    D2: ", solver.D2)
        println(io, "    D3: ", solver.D3)
    end
end

grid(solver::Solver) = grid(solver.D1)

@inline eachnode(solver::Solver) = Base.OneTo(length(grid(solver)))
@inline Base.real(::Solver{RealT}) where {RealT} = RealT

# Adapted from Trixi.jl
# https://github.com/trixi-framework/Trixi.jl/blob/75d8c67629562efd24b2a04e46d22b0a1f4f572c/src/solvers/dg.jl#L539
@inline function get_node_vars(q, equations, indices...)
    # There is a cut-off at `n == 10` inside of the method
    # `ntuple(f::F, n::Integer) where F` in Base at ntuple.jl:17
    # in Julia `v1.5`, leading to type instabilities if
    # more than ten variables are used. That's why we use
    # `Val(...)` below.
    # We use `@inline` to make sure that the `getindex` calls are
    # really inlined, which might be the default choice of the Julia
    # compiler for standard `Array`s but not necessarily for more
    # advanced array types such as `PtrArray`s, cf.
    # https://github.com/JuliaSIMD/VectorizationBase.jl/issues/55
    SVector(ntuple(@inline(v->q.x[v][indices...]), Val(nvariables(equations))))
end

@inline function set_node_vars!(q, q_node, equations, indices...)
    for v in eachvariable(equations)
        q.x[v][indices...] = q_node[v]
    end
    return nothing
end

function allocate_coefficients(mesh, equations, solver)
    return ArrayPartition(ntuple(_ -> zeros(real(solver), nnodes(mesh)),
                                 Val(nvariables(equations))))
end

function compute_coefficients!(q, func, t, mesh,
                               is_hyperbolic_appproximation::Val{false},
                               equations, solver)
    x = grid(solver)
    for i in eachnode(solver)
        q_node = func(x[i], t, equations, mesh)
        set_node_vars!(q, q_node, equations, i)
    end
end

# `set_approximation_variables!` is defined for equations, which are hyperbolic
# approximations, but below it is called for `equations` for which the compiler
# doesn't know from type annotations that these are hyperbolic approximations.
# Therefore, JET.jl fails. We provide a fallback definition here that does nothing.
set_approximation_variables!(q, mesh, equations, solver) = nothing

# For a hyperbolic approximation, we allow returning either the full set
# of variables or a reduced number determining only the limit system.
# In the second case, we compute the remaining variables using the default
# initialization specific to the equations.
function compute_coefficients!(q, func, t, mesh,
                               is_hyperbolic_appproximation::Val{true},
                               equations, solver)
    x = grid(solver)
    q_node = func(x[begin], t, equations, mesh)
    if length(q_node) == nvariables(equations)
        # full set of variables, proceed as usual
        compute_coefficients!(q, func, t, mesh, Val(false), equations, solver)
    else
        # reduced set of variables, fill remaining ones with NaN and
        # initialize them afterwards
        for i in eachnode(solver)
            q_node = func(x[i], t, equations, mesh)
            z = ntuple(_ -> convert(eltype(q_node), NaN),
                       Val(nvariables(equations) - length(q_node)))
            q_node = SVector(q_node..., z...)
            set_node_vars!(q, q_node, equations, i)
        end
        set_approximation_variables!(q, mesh, equations, solver)
    end
end

function calc_error_norms(q, t, initial_condition, mesh, equations,
                          solver)
    q_exact = similar(q)
    compute_coefficients!(q_exact, initial_condition, t, mesh,
                          is_hyperbolic_appproximation(equations), equations,
                          solver)
    l2_error = zeros(real(solver), nvariables(equations))
    linf_error = similar(l2_error)
    for v in eachvariable(equations)
        diff = q.x[v] - q_exact.x[v]
        l2_error[v] = integrate(abs2, diff, solver.D1) |> sqrt
        linf_error[v] = maximum(abs.(diff))
    end
    return l2_error, linf_error
end

function calc_sources!(dq, q, t, source_terms::Nothing, equations, solver)
    return nothing
end

function calc_sources!(dq, q, t, source_terms, equations, solver)
    x = grid(solver)
    for i in eachnode(solver)
        local_source = source_terms(get_node_vars(q, equations, i), x[i], t, equations)
        for v in eachvariable(equations)
            dq.x[v][i] += local_source[v]
        end
    end
    return nothing
end
