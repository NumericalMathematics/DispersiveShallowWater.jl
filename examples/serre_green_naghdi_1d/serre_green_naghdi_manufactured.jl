using OrdinaryDiffEqTsit5
using DispersiveShallowWater
using SummationByPartsOperators: upwind_operators, periodic_derivative_operator

###############################################################################
# Semidiscretization of the Serre-Green-Naghdi equations

# bathymetry_flat, bathymetry_mild_slope or bathymetry_variable
# be aware that for variable and mild_slope, at t=1.0 eta is constant in x
equations = SerreGreenNaghdiEquations1D(bathymetry_type = bathymetry_variable,
                                        gravity = 9.81)

initial_condition = initial_condition_manufactured
source_terms = source_terms_manufactured
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = 0.0
coordinates_max = 1.0
N = 256
mesh = Mesh1D(coordinates_min, coordinates_max, N)

#= for testing upwind operators use:
D1 = upwind_operators(periodic_derivative_operator;
                      derivative_order = 1, accuracy_order = 4,
                      xmin = xmin(mesh), xmax = xmax(mesh),
                      N = nnodes(mesh))
solver = Solver(D1) 
=#

# for testing central operators use:
acc_order = 4
solver = Solver(mesh, acc_order)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions,
                          source_terms = source_terms)

###############################################################################
# Create `ODEProblem` and run the simulation

tspan = (0.0, 1.2)
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 entropy_modified))

callbacks = CallbackSet(analysis_callback, summary_callback)
saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), abstol = 1e-12, reltol = 1e-12,
            save_everystep = false, callback = callbacks, saveat = saveat)
