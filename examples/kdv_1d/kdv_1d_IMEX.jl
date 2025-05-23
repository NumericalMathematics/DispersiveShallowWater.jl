using OrdinaryDiffEqTsit5
using DispersiveShallowWater
using SummationByPartsOperators: upwind_operators, periodic_derivative_operator

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
                          boundary_conditions = boundary_conditions)

tspan = (0.0, 5.0)
ode = semidiscretize(semi, tspan, no_splitform = false)

summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 waterheight))
callbacks = CallbackSet(analysis_callback, summary_callback)
saveat = range(tspan..., length = 100)

# alg = Rodas5() # not working because of https://github.com/SciML/OrdinaryDiffEq.jl/issues/2719
# alg = KenCarp4() # would need to add OrdinaryDiffEqSDIRK - can to that, but I dont want to bloat DSW.jl with package s
alg = Tsit5()
sol = solve(ode, alg, abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
