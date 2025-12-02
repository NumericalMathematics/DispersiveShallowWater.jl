using OrdinaryDiffEqRosenbrock
using DispersiveShallowWater
using SummationByPartsOperators: upwind_operators, periodic_derivative_operator

###############################################################################
# Semidiscretization of the KdV equation 

equations = KdVEquation1D(gravity = 9.81, D = 1.0)

initial_condition = initial_condition_manufactured
source_terms = source_terms_manufactured
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -1.0
coordinates_max = 1.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# Create solver with periodic SBP operators of accuracy order 3,
# which results in a 4th order accurate semi discretizations.
# We can set the accuracy order of the upwind operators to 3 since
# we only use central versions/combinations of the upwind operators.
D1_upwind = upwind_operators(periodic_derivative_operator;
                             derivative_order = 1, accuracy_order = 3,
                             xmin = xmin(mesh), xmax = xmax(mesh),
                             N = nnodes(mesh))
solver = Solver(D1_upwind)

semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions,
                          source_terms = source_terms)

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 waterheight, entropy))
callbacks = CallbackSet(analysis_callback, summary_callback)
saveat = range(tspan..., length = 100)

# The problem is very stiff. `Tsit5()` will cause a `maxiter` exceeded error.
alg = Rodas5()
sol = solve(ode, alg, abstol = 1e-12, reltol = 1e-12,
            save_everystep = false, callback = callbacks, saveat = saveat)
