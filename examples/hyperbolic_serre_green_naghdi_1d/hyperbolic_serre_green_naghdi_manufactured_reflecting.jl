using OrdinaryDiffEqTsit5
using DispersiveShallowWater

###############################################################################
# Semidiscretization of the hyperbolic Serre-Green-Naghdi equations

equations = HyperbolicSerreGreenNaghdiEquations1D(bathymetry_type = bathymetry_mild_slope,
                                                  lambda = 1.0e4,
                                                  gravity = 9.81)

initial_condition = initial_condition_manufactured_reflecting
source_terms = source_terms_manufactured_reflecting
boundary_conditions = boundary_condition_reflecting

# create homogeneous mesh
coordinates_min = 0.0
coordinates_max = 1.0
N = 128
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# create solver with SBP operators of accuracy order 4
accuracy_order = 4
D1 = derivative_operator(MattssonNordström2004();
                         derivative_order = 1, accuracy_order,
                         xmin = xmin(mesh), xmax = xmax(mesh), N = N)

solver = Solver(D1)
# semidiscretization holds all the necessary data structures for the spatial discretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions,
                          source_terms = source_terms)

###############################################################################
# Create `ODEProblem` and run the simulation
tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)
summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 1000,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 momentum, entropy))
callbacks = CallbackSet(analysis_callback, summary_callback)
saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), abstol = 1e-9, reltol = 1e-9,
            saveat = saveat, callback = callbacks)
plot(semi => sol, plot_initial = true, legend = false)


"""
accuracy_order = 2
####################################################################################################
l2
η                        v                        D                        w                        H
N    error     EOC       N    error     EOC       N    error     EOC       N    error     EOC       N    error     EOC
32   1.96e-02  -         32   6.18e-03  -         32   0.00e+00  -         32   1.11e+00  -         32   1.97e-02  -
64   6.60e-03  1.57      64   1.47e-03  2.07      64   0.00e+00  NaN       64   1.12e+00  -0.01     64   6.66e-03  1.57
128  2.27e-03  1.54      128  3.59e-04  2.03      128  0.00e+00  NaN       128  1.13e+00  -0.00     128  2.30e-03  1.53
256  7.93e-04  1.52      256  8.89e-05  2.02      256  0.00e+00  NaN       256  1.13e+00  -0.00     256  8.04e-04  1.52
512  2.79e-04  1.51      512  2.21e-05  2.01      512  0.00e+00  NaN       512  1.13e+00  -0.00     512  2.83e-04  1.51

mean           1.53      mean           2.03      mean           NaN       mean           -0.00     mean           1.53
----------------------------------------------------------------------------------------------------
linf
η                        v                        D                        w                        H
N    error     EOC       N    error     EOC       N    error     EOC       N    error     EOC       N    error     EOC
32   1.28e-01  -         32   9.91e-03  -         32   0.00e+00  -         32   1.72e+00  -         32   1.30e-01  -
64   6.46e-02  0.98      64   2.35e-03  2.08      64   0.00e+00  NaN       64   1.73e+00  -0.01     64   6.57e-02  0.99
128  3.25e-02  0.99      128  5.72e-04  2.04      128  0.00e+00  NaN       128  1.74e+00  -0.00     128  3.30e-02  0.99
256  1.63e-02  1.00      256  1.41e-04  2.02      256  0.00e+00  NaN       256  1.74e+00  -0.00     256  1.65e-02  1.00
512  8.15e-03  1.00      512  4.08e-05  1.79      512  0.00e+00  NaN       512  1.74e+00  -0.00     512  8.28e-03  1.00

mean           0.99      mean           1.98      mean           NaN       mean           -0.00     mean           0.99
----------------------------------------------------------------------------------------------------

accuracy_order = 4
####################################################################################################
l2
η                        v                        D                        w                        H
N    error     EOC       N    error     EOC       N    error     EOC       N    error     EOC       N    error     EOC
32   1.22e-02  -         32   5.39e-03  -         32   0.00e+00  -         32   1.14e+00  -         32   1.23e-02  -
64   2.53e-03  2.26      64   7.51e-04  2.84      64   0.00e+00  NaN       64   1.13e+00  0.01      64   2.57e-03  2.26
128  4.80e-04  2.40      128  9.78e-05  2.94      128  0.00e+00  NaN       128  1.13e+00  0.00      128  4.87e-04  2.40
256  8.77e-05  2.45      256  1.25e-05  2.97      256  0.00e+00  NaN       256  1.13e+00  0.00      256  8.89e-05  2.45
512  1.58e-05  2.47      512  2.09e-06  2.58      512  0.00e+00  NaN       512  1.13e+00  0.00      512  1.60e-05  2.47

mean           2.40      mean           2.83      mean           NaN       mean           0.00      mean           2.40
----------------------------------------------------------------------------------------------------
linf
η                        v                        D                        w                        H
N    error     EOC       N    error     EOC       N    error     EOC       N    error     EOC       N    error     EOC
32   9.73e-02  -         32   8.40e-03  -         32   0.00e+00  -         32   1.74e+00  -         32   9.91e-02  -
64   2.93e-02  1.73      64   1.21e-03  2.80      64   0.00e+00  NaN       64   1.74e+00  0.00      64   2.98e-02  1.73
128  7.94e-03  1.88      128  1.60e-04  2.91      128  0.00e+00  NaN       128  1.74e+00  0.00      128  8.06e-03  1.89
256  2.06e-03  1.95      256  2.77e-05  2.54      256  0.00e+00  NaN       256  1.74e+00  -0.00     256  2.09e-03  1.95
512  5.28e-04  1.96      512  3.78e-05  -0.45     512  0.00e+00  NaN       512  1.74e+00  0.00      512  5.32e-04  1.97

mean           1.88      mean           1.95      mean           NaN       mean           0.00      mean           1.89
----------------------------------------------------------------------------------------------------
"""