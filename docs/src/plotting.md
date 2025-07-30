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

DispersiveShallowWater.jl provides flexible plotting capabilities through [Plots.jl](https://github.com/JuliaPlots/Plots.jl) recipes. The plotting system supports various conversion functions, visualization options, and analysis tools. 

[Markie.jl](https://docs.makie.org/stable/) is not supported yet. [Contributions are welcome](https://github.com/NumericalMathematics/DispersiveShallowWater.jl/issues/220).

## Variable Conversion and Visualization Options

The plotting system supports different variable conversions and visualization options. You can plot conservative variables, specific physical quantities, and control what additional information is displayed:

```@example plotting
using Plots
# default(grid=true, box=:on, dpi=100, titlefont=font(16), linewidth=3, gridlinewidth=2, markersize=4, markerstrokewidth=2, xtickfontsize=14, ytickfontsize=14, xguidefontsize=16, yguidefontsize=16, ztickfontsize=14, zguidefontsize=16, legendfontsize=14) # hide

# Plot different variable representations

t = 13.37 # plot solution at (roughly) t = 13.37s
step_idx = argmin(abs.(saveat .- t)) # get the closed point to 13.37
p1 = plot(semi => sol, conversion = prim2prim, plot_bathymetry = false, 
          suptitle = "Primitive Variables", step = step_idx)
p2 = plot(semi => sol, conversion = prim2cons, plot_bathymetry = false,
          suptitle = "Conservative Variables", step = step_idx)
p3 = plot(semi => sol, conversion = waterheight_total, plot_bathymetry = true,
          suptitle = "Total Water Height", step = step_idx)
p4 = plot(semi => sol, conversion = velocity, plot_initial = true,
          suptitle = "Velocity with Initial Condition", step = step_idx)

plot(p1, p2, p3, p4, layout = (2, 2), size = (1000, 700))
savefig("variable_conversions.png") # hide
nothing # hide
```

![variable conversions](variable_conversions.png)

## Time Series Analysis at Spatial Points

You can analyze the temporal evolution of the solution at specific spatial locations. This is particularly useful for understanding wave propagation and local dynamics:

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

Using the analysis callback results, we can visualize how conserved quantities evolve over time:

```@example plotting
# Plot conservation properties
conservation_data = integrals(analysis_callback2)

p1 = plot(tstops(analysis_callback2), conservation_data.waterheight_total,
          title = "Water Mass Conservation", xlabel = "t", ylabel = "∫η",
          label = "Total water height", color = :blue)

p2 = plot(tstops(analysis_callback2), conservation_data.velocity,
          title = "Momentum Conservation", xlabel = "t", ylabel = "∫v", 
          label = "Total velocity", color = :red)

p3 = plot(tstops(analysis_callback2), conservation_data.entropy .- conservation_data.entropy[1],
          title = "Energy Conservation\n(with relaxation)", xlabel = "t", ylabel = "ΔE",
          label = "Energy change", color = :green)

plot(p1, p2, p3, layout = (3, 1), size = (500, 900))
savefig("conservation_analysis.png") # hide
nothing # hide
```

![conservation analysis](conservation_analysis.png)

## Error Analysis Visualization

Compare the numerical solution with the initial condition evaluated at time t. Note that this represents the analytical solution only if the initial condition function describes an exact solution that varies with time. You can either plot the errors manually or use the built-in plotting recipe with `plot(analysis_callback, what = (:errors,))`:

```@example plotting
# Calculate and plot errors over time
error_data = errors(analysis_callback2)
time_points = tstops(analysis_callback2)

p1 = plot(time_points, error_data.l2_error[1, :], 
          title = "L² Error Evolution", xlabel = "t", ylabel = "L² error", 
          label = "η")
plot!(p1, time_points, error_data.l2_error[2, :], 
      label = "v")

p2 = plot(time_points, error_data.linf_error[1, :], 
          title = "L∞ Error Evolution", xlabel = "t", ylabel = "L∞ error",
          label = "η")
plot!(p2, time_points, error_data.linf_error[2, :], 
      label = "v")

plot(p1, p2, layout = (1, 2), size = (800, 400))
savefig("error_analysis.png") # hide
nothing # hide
```

![error analysis](error_analysis.png)

The plotting system supports all standard Plots.jl features like custom color schemes, annotations, and interactive backends. For more advanced plotting options, consult the [Plots.jl documentation](https://docs.juliaplots.org/).
