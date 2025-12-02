using OrdinaryDiffEqLowStorageRK
using DispersiveShallowWater
using SummationByPartsOperators: MattssonNordström2004, derivative_operator

###############################################################################
# Semidiscretization of the hyperbolic Serre-Green-Naghdi equations

equations = HyperbolicSerreGreenNaghdiEquations1D(bathymetry_type = bathymetry_flat,
                                                  lambda = 500.0,
                                                  gravity = 9.81)

initial_condition = initial_condition_soliton
boundary_conditions = boundary_condition_reflecting

# create homogeneous mesh
coordinates_min = -50.0
coordinates_max = 50.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with SBP operators of accuracy order 2
accuracy_order = 2
D1 = derivative_operator(MattssonNordström2004();
                         derivative_order = 1, accuracy_order,
                         xmin = xmin(mesh), xmax = xmax(mesh), N = N)

solver = Solver(D1)

# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

###############################################################################
# Create `ODEProblem` and run the simulation
periods = 1.0
tspan = (0.0, periods * (xmax(mesh) - xmin(mesh)) / sqrt(1.2 * gravity(equations)))
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 energy_total_modified,
                                                                 entropy_modified))
relaxation_callback = RelaxationCallback(invariant = entropy_modified)
# Always put relaxation_callback before analysis_callback to guarantee conservation of the invariant
callbacks = CallbackSet(relaxation_callback, analysis_callback, summary_callback)

saveat = range(tspan..., length = 100)
# optimized time integration methods like this one are much more efficient
# for stiff problems (λ big) than standard methods like Tsit5()
alg = RDPK3SpFSAL35()
sol = solve(ode, alg, abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
