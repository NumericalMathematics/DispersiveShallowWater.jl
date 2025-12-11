# Basic Example

## Introduction

In this tutorial we describe how to numerically solve the BBM-BBM (Benjamin-Bona-Mahony) equations with variable bottom topography in one dimension. The equations describe a dispersive shallow water model,
i.e. they extend the well-known shallow water equations in the sense that dispersion is modeled. The shallow water equations are a system of first
order hyperbolic partial differential equations that can be written in the form of a balance law. In contrast, the BBM-BBM equations additionally
include third-order mixed derivatives. In primitive variables ``q = (\eta, v)`` they can be written as:

```math
\begin{aligned}
  \eta_t + ((\eta + D)v)_x - \frac{1}{6}(D^2\eta_{xt})_x &= 0,\\
  v_t + g\eta_x + \left(\frac{1}{2}v^2\right)_x - \frac{1}{6}(D^2v_t)_{xx} &= 0.
\end{aligned}
```

All equations in [DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl) follow the same variable conventions: ``\eta = h + b`` describes the total water height, ``h`` the water height above the bottom topography (bathymetry), ``b = \eta_0 - D`` the bathymetry, and ``v`` the velocity in horizontal direction. The reference water height ``\eta_0`` (also called still water height) and gravitational acceleration ``g`` are used consistently across all equations. For the BBM-BBM equations specifically, ``\eta_0`` is typically set to 0.

A sketch of the water height and bathymetry can be found below.

![water height and bathymetry](bathymetry.png)

## [Getting started](@id basic_example)

In order to conduct a numerical simulation with [DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl), we perform the following steps.

First, we load the necessary libraries:

```@example overview
using DispersiveShallowWater, OrdinaryDiffEqTsit5
```

## Define physical setup

As a first step of a numerical simulation, we define the physical setup we want to solve. This includes the set of equations, potentially
including physical parameters, initial and boundary conditions as well as the domain. In the following example, the initial condition
describes a traveling wave that moves towards a beach, which is modeled by a linearly increasing bathymetry.

```@example overview
equations = BBMBBMEquations1D(bathymetry_type = bathymetry_variable, gravity = 9.81)

function initial_condition_shoaling(x, t, equations::BBMBBMEquations1D, mesh)
    A = 0.07 # amplitude of wave
    x0 = -30 # initial center
    eta = A * exp(-0.1*(x - x0)^2)
    v = 0
    D = x <= 0.0 ? 0.7 : 0.7 - 1/50 * x
    return SVector(eta, v, D)
end

initial_condition = initial_condition_shoaling
boundary_conditions = boundary_condition_periodic

coordinates_min = -130.0
coordinates_max = 20.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)
```

The first line specifies that we want to solve the BBM-BBM equations with variable bathymetry using a gravitational acceleration of ``g = 9.81``.
Afterwards, we define the initial condition, which is described as a function with the spatial variable `x`, the time `t`, the `equations`, and
a `mesh` as parameters.

If an analytical solution is available, the time variable `t` can be used, and the initial condition can serve as an analytical solution to be compared with the numerical solution. Otherwise, you can just keep the time variable unused.

An initial condition in [DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl) is supposed to return an `SVector` holding the values for each of the unknown variables. Since the bathymetry is
treated as a variable (with time derivative 0) for convenience, we need to provide the value for the primitive variables `eta` and `v` as well as for `D`.

Next, we choose periodic boundary conditions. [DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl) also supports reflecting boundary conditions for some but not all equations. For more information see the [Dispersive Shallow Water Models overview](@ref eq_overview).

Lastly, we define the physical domain as the interval from -130 to 20 and we choose 512 nodes. The node distribution of the mesh depends on the chosen SBP operators. For classical finite difference operators this is homogeneous, i.e. the distance between consecutive nodes is constant. We choose the left boundary very far to the left in order to avoid interactions between left- and right-traveling waves. This prevents unwanted wave interference that could occur when waves wrap around due to the periodic boundary conditions.

## Define numerical solver

In the next step, we build a [`Semidiscretization`](@ref) that bundles all ingredients for the spatial discretization of the model. Especially, we need to define a [`Solver`](@ref).

The simplest way to define a solver when working with [`boundary_condition_periodic`](@ref) is to call the constructor by providing the mesh and a desired order of accuracy.

In the following example, we use an accuracy order of 4. The default constructor simply creates periodic first-, second-, and third-derivative central finite difference summation-by-parts (SBP) operators of the provided order of accuracy.

How to use other summation-by-parts operators, is described in the section on [how to customize the solver](@ref customize_solver). Note that for non-periodic boundary conditions, the solver also needs to be created with non-periodic
operators, see, e.g. [examples/bbm\_bbm\_1d/bbm\_bbm\_1d\_basic\_reflecting.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl/blob/main/examples/bbm_bbm_1d/bbm_bbm_1d_basic_reflecting.jl).

```@example overview
accuracy_order = 4
solver = Solver(mesh, accuracy_order)

semi = Semidiscretization(mesh, equations, initial_condition, solver, boundary_conditions = boundary_conditions)
```

Finally, we put the `mesh`, the `equations`, the `initial_condition`, the `solver`, and the `boundary_conditions` together in a semidiscretization `semi`.

## Solve system of ordinary differential equations

Once we have obtained a semidiscretization, we can solve the resulting system of ordinary differential equations. To do so, we specify the time interval that
we want to simulate and obtain an `ODEProblem` from the [SciML ecosystem for ordinary differential equations](https://diffeq.sciml.ai/latest/) by calling
[`semidiscretize`](@ref) on the semidiscretization and the time span.

Finally, the `ode` can be `solve`d using the interface from [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl). This means we can specify a time-stepping
scheme we want to use, the tolerances for the adaptive time-stepping, and the time values, where the solution values should be saved. In this case, we use the adaptive
explicit Runge-Kutta method `Tsit5` by Tsitouras of order 5(4), which is implemented in the subpackage OrdinaryDiffEqTsit5.jl. If you want to use other time-stepping
schemes, you can install the respective subpackage or the whole package OrdinaryDiffEq.jl, which will install every available solver.
Here, we save the solution at 100 equidistant points in time.

```@example overview
tspan = (0.0, 25.0)
ode = semidiscretize(semi, tspan)

saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7, saveat = saveat)
```

After solving the equations, `sol` contains the solution for each of the—in this case—three variables at every spatial point for each of the 100 points in time.

## [Visualize results](@id visualize_results)

After running the simulation, the results can be visualized using [Plots.jl](https://github.com/JuliaPlots/Plots.jl), which needs to be imported first. Then, we can
plot the solution at the final time by calling `plot` on a `Pair` of the `Semidiscretization` and the corresponding `ODESolution` `sol`. The result is depicted in the following picture.

```@example overview
using Plots, Printf
# default(grid=true, box=:on, dpi=100, titlefont=font(16), linewidth=3, gridlinewidth=2, markersize=4, markerstrokewidth=2, xtickfontsize=14, ytickfontsize=14, xguidefontsize=16, yguidefontsize=16, ztickfontsize=14, zguidefontsize=16, legendfontsize=14) # hide
plot(semi => sol, plot_title = @sprintf "BBM-BBM equations at t = %.2f" last(tspan))
savefig("shoaling_solution.png") # hide
nothing # hide
```

![shoaling solution](shoaling_solution.png)

By default, this will plot the bathymetry, but not the initial (analytical) solution.

You can adjust this by passing the boolean values `plot_bathymetry` (if `true`, always plot bathymetry in the first subplot) and `plot_initial`. Note that `plot_initial = true` will evaluate and plot the initial condition function at the same time `t` as the numerical solution being displayed (the final time by default). This means if your initial condition function represents an analytical solution, setting `plot_initial = true` will plot the analytical solution at that specific time for comparison.

Plotting an animation over time can, e.g., be done by the following command, which uses `step` to plot the solution at a specific time step. Here `conversion = waterheight_total` makes it so that we only look at the total water height ``\eta`` and not also the velocity ``v``. More on tutorials for plotting can be found in the chapter [Plotting Simulation Results](@ref plotting).

```@example overview
anim = @animate for step in 1:length(sol.u)
    plot(semi => sol, plot_initial = true, conversion = waterheight_total, step = step, xlims = (-50, 20), ylims = (-0.8, 0.1),
         plot_title = @sprintf "BBM-BBM equations at t = %.2f" sol.t[step])
end
gif(anim, "shoaling_solution.gif", fps = 25)
nothing # hide
```

![shoaling solution](shoaling_solution.gif)

It is also possible to plot the solution variables at a fixed spatial point over time by calling `plot(semi => sol, x)` for some `x`-value, see the chapter [Plotting Simulation Results](@ref plotting) in this documentation for some examples.

## More examples

More examples sorted by the simulated equations can be found in the [examples/](https://github.com/NumericalMathematics/DispersiveShallowWater.jl/tree/main/examples) subdirectory.

## [Plain program](@id overview-plain-program)

Here follows a version of the program without any comments.

```Julia
using DispersiveShallowWater, OrdinaryDiffEqTsit5

equations = BBMBBMEquations1D(bathymetry_type = bathymetry_variable, gravity = 9.81)

function initial_condition_shoaling(x, t, equations::BBMBBMEquations1D, mesh)
    A = 0.07 # amplitude of wave
    x0 = -30 # initial center
    eta = A * exp(-0.1*(x - x0)^2)
    v = 0
    D = x <= 0.0 ? 0.7 : 0.7 - 1/50 * x
    return SVector(eta, v, D)
end

initial_condition = initial_condition_shoaling
boundary_conditions = boundary_condition_periodic

coordinates_min = -130.0
coordinates_max = 20.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

solver = Solver(mesh, 4)

semi = Semidiscretization(mesh, equations, initial_condition, solver, boundary_conditions = boundary_conditions)

tspan = (0.0, 25.0)
ode = semidiscretize(semi, tspan)

saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7, saveat = saveat)

using Plots
plot(semi => sol)

anim = @animate for step in 1:length(sol.u)
    plot(semi => sol, plot_initial = true, conversion = waterheight_total, step = step, xlims = (-50, 20), ylims = (-0.8, 0.1))
end
gif(anim, "shoaling_solution.gif", fps = 25)
```
