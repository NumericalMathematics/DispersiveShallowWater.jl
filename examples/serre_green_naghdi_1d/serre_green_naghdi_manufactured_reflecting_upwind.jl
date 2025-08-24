using OrdinaryDiffEqTsit5
using DispersiveShallowWater
using SummationByPartsOperators: Mattsson2017, upwind_operators

###############################################################################
# Semidiscretization of the Serre-Green-Naghdi equations

equations = SerreGreenNaghdiEquations1D(bathymetry_type = bathymetry_mild_slope,
                                        gravity = 9.81)

initial_condition = initial_condition_manufactured_reflecting
source_terms = source_terms_manufactured_reflecting
boundary_conditions = boundary_condition_reflecting

# create homogeneous mesh
coordinates_min = 0.0
coordinates_max = 1.0
N = 64
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with SBP operators of accuracy order 2
accuracy_order = 2
D1 = upwind_operators(Mattsson2017, derivative_order = 1, accuracy_order = accuracy_order,
                      xmin = xmin(mesh), xmax = xmax(mesh),
                      N = nnodes(mesh))
solver = Solver(D1)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver;
                          boundary_conditions, source_terms)

###############################################################################
# Create `ODEProblem` and run the simulation
tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_integrals = (waterheight_total,
                                                                 velocity, entropy_modified,
                                                                 entropy))
callbacks = CallbackSet(analysis_callback, summary_callback)
sol = solve(ode, Tsit5(); abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks)
