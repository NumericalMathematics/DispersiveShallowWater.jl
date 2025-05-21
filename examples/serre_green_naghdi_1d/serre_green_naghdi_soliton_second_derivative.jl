using OrdinaryDiffEqTsit5
using DispersiveShallowWater
using SummationByPartsOperators: Mattsson2012, derivative_operator,
                                 var_coef_derivative_operator

###############################################################################
# Semidiscretization of the Serre-Green-Naghdi equations

equations = SerreGreenNaghdiEquations1D(bathymetry_type = bathymetry_flat,
                                        gravity = 9.81)

# initial_condition_convergence_test can only be used to get reasonable errors
# for periodic boundary conditions - but we can nevertheless compute the
# evolution of the soliton with reflecting boundary conditions
initial_condition = initial_condition_convergence_test
boundary_conditions = boundary_condition_reflecting

# create homogeneous mesh
coordinates_min = -50.0
coordinates_max = 50.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with SBP operators of accuracy order 4
accuracy_order = 4
D1 = derivative_operator(Mattsson2012();
                         derivative_order = 1, accuracy_order,
                         xmin = xmin(mesh), xmax = xmax(mesh), N = N)
D2 = var_coef_derivative_operator(Mattsson2012(),
                                  2, accuracy_order,
                                  xmin(mesh), xmax(mesh), N, one)
solver = Solver(D1, D2)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation
# FIXME
factor = 0.5
tspan = (0.0, factor * (xmax(mesh) - xmin(mesh)) / sqrt(1.2 * equations.gravity)) # one period
# tspan = (0.0, (xmax(mesh) - xmin(mesh)) / sqrt(1.2 * equations.gravity)) # one period
ode = semidiscretize(semi, tspan)
ode.u0.x[2][begin] = 0;
ode.u0.x[2][end] = 0; # FIXME
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 entropy_modified))
# callbacks = CallbackSet(analysis_callback, summary_callback)
# FIXME: The current way to define a pointwise modified energy does not work
#        immediately since we define the derivative term using a quadratic form in
#        the velocity
callbacks = CallbackSet(summary_callback)

# saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(); abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks)
