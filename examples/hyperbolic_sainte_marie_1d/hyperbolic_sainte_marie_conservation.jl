# This elixir contains an artificial setup that can be used to check the
# conservation properties of the equations and numerical methods as well as
# a possible directional bias (if the velocity is set to zero). See
# - Hendrik Ranocha and Mario Ricchiuto (2024)
#   Structure-preserving approximations of the Serre-Green-Naghdi
#   equations in standard and hyperbolic form
#   [arXiv: 2408.02665](https://arxiv.org/abs/2408.02665)

using OrdinaryDiffEqLowStorageRK
using DispersiveShallowWater

#=
FIXME: Test this using
julia> using Revise; using DispersiveShallowWater, Plots, OrdinaryDiffEqVerner

julia> trixi_include("examples/hyperbolic_sainte_marie_1d/hyperbolic_sainte_marie_conservation.jl", tol = 1.0e-12, alg = Vern9());

julia> integrals(analysis_callback).waterheight_total |> x -> (x[end] - x[1]) / x[1]
0.0

julia> integrals(analysis_callback).momentum |> x -> (x[end] - x[1]) / x[1]
2.9432057444810533e-15

julia> integrals(analysis_callback).entropy_modified |> x -> (x[end] - x[1]) / x[1]
-1.064590038909088e-15
=#

###############################################################################
# Semidiscretization of the hyperbolic Sainte-Marie equations

equations = HyperbolicSainteMarieEquations1D(bathymetry_type = bathymetry_flat,
                                             gravity = 9.81,
                                             eta0 = 1.0,
                                             alpha = 3.0)

function initial_condition_conservation_test(x, t,
                                             equations::HyperbolicSainteMarieEquations1D,
                                             mesh)
    eta = 1 + exp(-x^2)
    v = 1.0e-2 # set this to zero to test a directional bias
    b = 0.0 # FIXME: 0.25 * cospi(x / 75)

    # We use the feature that we can only return the physical variables
    # used by the `hyperbolic_approximation_limit`, i.e., the
    # `SainteMarieEquations1D.`
    D = equations.eta0 - b
    return SVector(eta, v, D)
end

# create homogeneous mesh
coordinates_min = -150.0
coordinates_max = +150.0
N = 1_000
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with periodic SBP operators of accuracy order 2
accuracy_order = 2
solver = Solver(mesh, accuracy_order)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations,
                          initial_condition_conservation_test, solver;
                          boundary_conditions = boundary_condition_periodic)

###############################################################################
# Create `ODEProblem` and run the simulation
tspan = (0.0, 35.0)
ode = semidiscretize(semi, tspan)

# The callbacks support an additional `io` argument to write output to a file
# or any other IO stream. The default is stdout. We use this here to enable
# setting it to `devnull` to benchmark the full simulation including the time
# to compute the errors etc. but without the time to write the output to the
# terminal.
io = stdout
summary_callback = SummaryCallback(io)
analysis_callback = AnalysisCallback(semi; interval = 50, io,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 momentum,
                                                                 entropy_modified))

callbacks = CallbackSet(analysis_callback, summary_callback)

# optimized time integration methods like this one are much more efficient
# for stiff problems (α big) than standard methods like Tsit5()
alg = RDPK3SpFSAL35()
tol = 1.0e-4
sol = solve(ode, alg;
            reltol = tol, abstol = tol,
            save_everystep = false, callback = callbacks)
