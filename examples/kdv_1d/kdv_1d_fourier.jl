using OrdinaryDiffEqTsit5
using DispersiveShallowWater
using SummationByPartsOperators: fourier_derivative_operator

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

# create solver with Fourier SBP operators
D1 = fourier_derivative_operator(xmin(mesh), xmax(mesh), nnodes(mesh))
D3 = D1^3

solver = Solver(D1, nothing, D3)

semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

tspan = (0.0, 5.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 waterheight))
callbacks = CallbackSet(analysis_callback, summary_callback)
saveat = range(tspan..., length = 100)

# Fourier SBP operators are not working with ForwardDiff.jl, so use an explicit method.
alg = Tsit5()
sol = solve(ode, alg, abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
