using OrdinaryDiffEqTsit5
using DispersiveShallowWater

###############################################################################
# Semidiscretization of the KdV equation 

equations = KdVEquation1D(gravity = 9.81, D = 1.0)
initial_condition = initial_condition_convergence_test
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -50.0
coordinates_max = 50.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with periodic SBP operators of accuracy order 4
accuracy_order = 4
solver = Solver(mesh, accuracy_order)

semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

tspan = (0.0, 5.0)
ode = semidiscretize(semi, tspan, split_ode = Val{false}()) # no IMEX for now

summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 waterheight))
callbacks = CallbackSet(analysis_callback, summary_callback)
saveat = range(tspan..., length = 100)

sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
