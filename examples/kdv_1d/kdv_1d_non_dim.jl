using OrdinaryDiffEqTsit5
using DispersiveShallowWater
using SummationByPartsOperators: upwind_operators, periodic_derivative_operator

###############################################################################
# Semidiscretization of the KdV equation 

function initial_condition_non_dimensional(x, t, equations::KdVEquation1D, mesh)
    c = 1 / 3
    x_t = mod(x - c * t - xmin(mesh), xmax(mesh) - xmin(mesh)) + xmin(mesh)
    u = 3 * c * sech(sqrt(9 * c) / 6 * x_t)^2
    return SVector(u)
end

function initial_condition_non_dimensional_converted(x, t, equations::KdVEquation1D, mesh)
    u = initial_condition_non_dimensional(x, t, equations, mesh)
    return nondim2prim(u, equations) # return eta
end

# Parameters g = 4/27 and D = 3.0 are needed for conversion to non-dimensional variables
equations = KdVEquation1D(gravity = 4 / 27, D = 3.0)
initial_condition = initial_condition_non_dimensional_converted
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

tspan = (0.0, 50.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 waterheight, entropy))
callbacks = CallbackSet(analysis_callback, summary_callback)
saveat = range(tspan..., length = 100)

sol = solve(ode, Tsit5(), abstol = 1e-8, reltol = 1e-8,
            save_everystep = false, callback = callbacks, saveat = saveat)

# Plot the solution transformed back to non dimensional variables.
# plot(semi => sol, conversion = prim2nondim, plot_initial = true)
