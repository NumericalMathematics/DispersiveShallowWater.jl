```@example plotting
using DispersiveShallowWater, OrdinaryDiffEqTsit5# hide
equations = BBMBBMEquations1D(gravity = 9.81)# hide
# hide
initial_condition = initial_condition_convergence_test# hide
# hide
# Create mesh and solver# hide
coordinates_min = -10.0# hide
coordinates_max = 10.0# hide
N = 128# hide
mesh = Mesh1D(coordinates_min, coordinates_max, N)# hide
solver = Solver(mesh, 4)# hide
# hide
# Create semidiscretization# hide
semi = Semidiscretization(mesh, equations, initial_condition, solver,# hide
                          boundary_conditions = boundary_condition_periodic)# hide
                          tspan = (0.0, 20.0)# hide
ode = semidiscretize(semi, tspan)# hide
# hide
analysis_callback = AnalysisCallback(semi; interval = 10,# hide
                                     extra_analysis_errors = (:conservation_error,),# hide
                                     extra_analysis_integrals = (waterheight_total,# hide
                                                                 velocity, entropy),# hide
                                     io = devnull)# hide
# hide
saveat = range(tspan..., length = 100)# hide
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,# hide
            save_everystep = false, callback = analysis_callback, saveat = saveat)# hide
analysis_callback2 = AnalysisCallback(semi; interval = 10,# hide
                                     extra_analysis_errors = (:conservation_error,),# hide
                                     extra_analysis_integrals = (waterheight_total,# hide
                                                                 velocity, entropy),# hide
                                     io = devnull)# hide
relaxation_callback = RelaxationCallback(invariant = entropy)# hide
callbacks = CallbackSet(relaxation_callback, analysis_callback2)# hide
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,# hide
            save_everystep = false, callback = callbacks, saveat = saveat)# hide

nothing # hide
```

# [Plotting Simulation Results](@id plotting)

[DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl) provides flexible plotting capabilities through [Plots.jl](https://github.com/JuliaPlots/Plots.jl) recipes. The plotting system supports various conversion functions, visualization options, and analysis tools.

[Makie.jl](https://docs.makie.org/stable/) is not supported yet. [Contributions are welcome](https://github.com/NumericalMathematics/DispersiveShallowWater.jl/issues/220).

## Variable Conversion and Visualization Options

The plotting system supports different variable conversions and visualization options. By default, the physical variables (returned by [`prim2phys`](@ref)) are plotted. For hyperbolic approximations like [`HyperbolicSerreGreenNaghdiEquations1D`](@ref), this means auxiliary variables are not plotted by default – only the physical variables from the limit system.

You can plot conservative variables, specific physical quantities, and control what additional information is displayed:

```@example plotting
using Plots
# default(grid=true, box=:on, dpi=100, titlefont=font(16), linewidth=3, gridlinewidth=2, markersize=4, markerstrokewidth=2, xtickfontsize=14, ytickfontsize=14, xguidefontsize=16, yguidefontsize=16, ztickfontsize=14, zguidefontsize=16, legendfontsize=14) # hide

# Plot different variable representations

t = 13.37 # plot solution at (roughly) t = 13.37s
step_idx = argmin(abs.(saveat .- t)) # get the closest point to 13.37
p1 = plot(semi => sol, conversion = prim2prim, plot_bathymetry = false, plot_analytical = true,
          suptitle = "Primitive Variables with Analytical Solution", step = step_idx)
p2 = plot(semi => sol, conversion = prim2cons, plot_bathymetry = false,
          suptitle = "Conservative Variables", step = step_idx)
p3 = plot(semi => sol, conversion = waterheight_total, plot_bathymetry = true,
          suptitle = "Total Water Height", step = step_idx)
p4 = plot(semi => sol, conversion = velocity, plot_initial = true, plot_bathymetry = false,
          suptitle = "Velocity with Initial Condition", step = step_idx)

plot(p1, p2, p3, p4, layout = (2, 2), size = (1000, 700))
savefig("variable_conversions.png") # hide
nothing # hide
```

![variable conversions](variable_conversions.png)

Note that the argument `plot_analytical = true` plots the initial condition evaluated at the selected time step, which means
that the analytical solution is plotted if the initial condition function describes an exact solution that varies with time.
Similarly, `plot_initial = true` always plots the initial condition at time `t = first(tspan)`.

## Time Series Analysis at Spatial Points

You can analyze the temporal evolution of the solution at specific spatial locations by passing a spatial point as second argument. This is particularly useful for understanding wave propagation and local dynamics:

```@example plotting
# Analyze solution at a single spatial point
x_location = 0.0

p1 = plot(semi => sol, x_location, conversion = waterheight_total,
          suptitle = "Water Height Evolution at x = $x_location")

p2 = plot(semi => sol, x_location, conversion = velocity,
          suptitle = "Velocity Evolution at x = $x_location")

plot(p1, p2, layout = (1, 2), size = (800, 400))
savefig("time_series_analysis.png") # hide
nothing # hide
```

![time series analysis](time_series_analysis.png)

## Energy and Momentum Evolution

Using the analysis callback results, we can visualize how conserved quantities evolve over time. The built-in plotting recipe shows the **change** of each invariant from its initial value:

```@example plotting
# Plot conservation properties using the built-in recipe
# This shows how each conserved quantity deviates from its initial value over time
plot(analysis_callback2, exclude = (:velocity,))
savefig("conservation_analysis.png") # hide
nothing # hide
```

![conservation analysis](conservation_analysis.png)

The plot shows the change of invariants over time (i.e., `integral(t) - integral(t=0)`). The `exclude` parameter allows you to hide specific quantities from the plot - here we exclude velocity to focus on other conserved quantities.

## Error Analysis Visualization

The error analysis compares the numerical solution with the initial condition evaluated at time t. Note that this represents the analytical solution only if the initial condition function describes an exact solution that varies with time:

```@example plotting
# Plot error evolution using the built-in recipe
# The 'what' parameter controls what gets plotted: (:integrals,), (:errors,), or both
# The 'exclude' parameter removes specific error types from the plot
plot(analysis_callback2, what = (:errors,), exclude = (:conservation_error,))
savefig("error_analysis.png") # hide
nothing # hide
```

![error analysis](error_analysis.png)

The `what = (:errors,)` parameter tells the plotting recipe to show errors instead of the default invariants. The errors plotted are the **total errors summed over all variables** (``L^2`` and ``L^\infty`` norms). The `exclude = (:conservation_error,)` parameter removes the conservation error from the plot, focusing only on the discretization errors (``L^2`` and ``L^\infty``).

The plotting system supports all standard [Plots.jl](https://github.com/JuliaPlots/Plots.jl) features like custom color schemes, annotations, and interactive backends. For more advanced plotting options, consult the [Plots.jl documentation](https://docs.juliaplots.org/).
