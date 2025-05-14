using OrdinaryDiffEqTsit5
using DispersiveShallowWater
using SummationByPartsOperators:  periodic_derivative_operator, upwind_operators

###############################################################################
# Semidiscretization of the KdV equation 

equations = KdVEquation1D(gravity = 9.81, D = 1.0)
initial_condition = initial_condition_convergence_test
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -50.0
coordinates_max = 50.0
N = 1024
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with periodic SBP operators of accuracy order 4
accuracy_order = 4

          

# using a wide stencil third derivative operator is significantly faster than using
# the upwind operator for the third derivative for a given N. However there are considerably more oscillations.
# This can be improved by increasing N, which in fact leads to it being slower again.
# TODO: This text maybe in the docs somewhere?
D1 = periodic_derivative_operator(1, accuracy_order, mesh.xmin, mesh.xmax, mesh.N)
D3 = periodic_derivative_operator(3, accuracy_order, mesh.xmin, mesh.xmax, mesh.N)

solver = Solver(D1, D3)




semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)                         

summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),)
callbacks = nothing
saveat = range(tspan..., length = 100)

@btime sol = solve(ode, Tsit5(), abstol = 1e-9, reltol = 1e-9,
            save_everystep = false, callback = callbacks, saveat = saveat)

plot(semi => sol)


"""
OG:
@btime = 507.327 ms (197787 allocations: 170.53 MiB)
With Cache:
@btime = 15.819 ms (3981 allocations: 542.41 KiB)
for Tsit5(), abstol = 1e-9, reltol = 1e-9
with 
(gravity = 9.81, D = 1.0)
initial_condition = initial_condition_convergence_test
boundary_conditions = boundary_condition_periodic
coordinates_min = -50.0
coordinates_max = 50.0
N = 512
"""