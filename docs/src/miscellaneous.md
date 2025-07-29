# Miscellaneous

# Callbacks

Callbacks provide additional functionality during simulations, such as monitoring solution properties, analyzing errors, or ensuring conservation of physical quantities. DispersiveShallowWater.jl implements three main callback types that can be used individually or in combination to enhance simulation analysis and performance monitoring.

When using multiple callbacks simultaneously, combine them using a `CallbackSet`:

```julia
callbacks = CallbackSet(analysis_callback, summary_callback)
sol = solve(ode, Tsit5(), callback = callbacks)
```

## Summary Callback

The [`SummaryCallback`](@ref) provides performance profiling information at the end of a simulation. It tracks computational time spent in different parts of the code and memory allocations, giving insights into the computational efficiency of the simulation.

The callback automatically prints a detailed timing breakdown showing:

- Total simulation time and memory allocations
- Time spent in different computational sections
- Number of function calls and average execution time per call

```julia
summary_callback = SummaryCallback()
sol = solve(ode, Tsit5(), callback = summary_callback)
```

At the end of the simulation, the callback will display output similar to:

```
───────────────────────────────────────────────────────────────────────────────────────────
              DispersiveSWE                       Time                    Allocations
                                         ───────────────────────   ────────────────────────
            Tot / % measured:                 440ms /  98.4%            867MiB /  99.9%

Section                          ncalls     time    %tot     avg     alloc    %tot      avg
───────────────────────────────────────────────────────────────────────────────────────────
rhs!                              1.62k    420ms   97.0%   258μs    866MiB  100.0%   546KiB
  solving elliptic system         1.62k    247ms   57.2%   152μs    357MiB   41.2%   225KiB
  assembling elliptic operator    1.62k    167ms   38.6%   103μs    509MiB   58.8%   321KiB
  hyperbolic terms                1.62k   3.49ms    0.8%  2.15μs     0.00B    0.0%    0.00B
  ~rhs!~                          1.62k   1.46ms    0.3%   902ns   1.55KiB    0.0%    0.98B
  source terms                    1.62k   42.0μs    0.0%  25.9ns     0.00B    0.0%    0.00B
analyze solution                      3   13.2ms    3.0%  4.39ms    147KiB    0.0%  49.1KiB
───────────────────────────────────────────────────────────────────────────────────────────
```


## Analysis Callback

The [`AnalysisCallback`](@ref) monitors solution quality and physical properties during the simulation. It computes error norms and tracks conservation of important physical quantities at specified time intervals.

### Setting up the Analysis Callback

First, let's set up a basic simulation using the BBM-BBM equations:

```@example callback
using DispersiveShallowWater, OrdinaryDiffEqTsit5

# Define the physical setup
equations = BBMBBMEquations1D(gravity = 9.81)

initial_condition = initial_condition_convergence_test

# Create mesh and solver
coordinates_min = -10.0
coordinates_max = 10.0
N = 128
mesh = Mesh1D(coordinates_min, coordinates_max, N)
solver = Solver(mesh, 4)

# Create semidiscretization
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_condition_periodic)
nothing # hide
```

### Error Analysis and Conservation Monitoring

The analysis callback computes ``L^2`` and ``L^\infty`` errors by comparing the numerical solution to the initial condition at time ``t`` (which can be the analytical solution, if available). Additional error types can be specified using the `extra_analysis_errors` parameter, and physical quantities can be monitored using `extra_analysis_integrals`.

The conservation error measures the temporal change of conserved quantities. For the BBM-BBM equations, important conserved quantities include the total water mass (integral of water height `h`), the total momentum (integral of `v` for flat bathymetry), and the [`entropy`](@ref). The specific form of the entropy varies between different equation systems. For the BBM-BBM equations, the entropy is:

```math
\mathcal E(t; \eta, v) = \frac{1}{2}\int_\Omega g\eta^2 + (\eta + D)v^2 \, dx
```

where ``\eta`` is the total water height and ``D`` is the still-water depth.

```@example callback
tspan = (0.0, 20.0)
ode = semidiscretize(semi, tspan)

analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 velocity, entropy),
                                     io = devnull)

saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = analysis_callback, saveat = saveat)
nothing # hide
```

The recorded errors and integrals can be accessed as `NamedTuple`s using [`errors(analysis_callback)`](@ref) and [`integrals(analysis_callback)`](@ref).

### Visualizing Conservation Properties

The temporal evolution of monitored quantities can be visualized by plotting the analysis callback:

```@example callback
using Plots
plot(analysis_callback)
savefig("analysis_callback.png") # hide
nothing # hide
```

![analysis callback](analysis_callback.png)

The plot shows that linear invariants such as the total water mass and total velocity are conserved exactly. However, nonlinear invariants such as the entropy may exhibit small growth over time. This occurs because standard time integration methods do not necessarily preserve nonlinear invariants, even when the spatial discretization is conservative.

## Relaxation Callback

To obtain entropy-conserving time-stepping schemes, DispersiveShallowWater.jl uses the relaxation method introduced in [^Ketcheson2019] and further developed in [^RanochaSayyariDalcinParsaniKetcheson2020]. The relaxation method is implemented as a [`RelaxationCallback`](@ref), which takes a function representing the conserved quantity as the keyword argument `invariant`. This callback modifies the time step to maintain conservation of a specified quantity.

### Entropy-Conserving Time Integration

To achieve exact conservation of the entropy, we add a relaxation callback to the simulation:

```@example callback
analysis_callback2 = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 velocity, entropy),
                                     io = devnull)
relaxation_callback = RelaxationCallback(invariant = entropy)

# Important: RelaxationCallback must come before AnalysisCallback
callbacks = CallbackSet(relaxation_callback, analysis_callback2)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)
nothing # hide
```

!!! note "Callback Ordering"
    When using both `RelaxationCallback` and `AnalysisCallback`, the relaxation callback must be placed first in the `CallbackSet`. This ensures that the analysis callback monitors the solution after the relaxation step has been applied.

The relaxation method modifies each time step by finding an optimal relaxation parameter that preserves the specified invariant exactly. This results in entropy conservation up to machine precision: NOT MACHINE PRECISION. WHY? Something strange with BBMBBM it feels like.

```@example callback
plot(analysis_callback2, ylims = (-2e-12, 2e-12))
savefig("analysis_callback_relaxation.png") # hide
nothing # hide
```

![analysis callback relaxation](analysis_callback_relaxation.png)

The plot demonstrates that with the relaxation callback, the entropy is conserved to machine precision throughout the simulation, providing a fully discrete structure-preserving numerical scheme.

## Relaxation

The relaxation method conserves nonlinear invariants up to machine precision in both the spatial and temporal discretization with minimal computational overhead. This can improve solution stability and accuracy, often reducing error growth over time from quadratic to linear.

The following comparison shows error growth with and without the relaxation method using the above simulation setup:

```@example callback
plot(errors(analysis_callback).l2_error[1, :], label = "without relaxation")
plot!(errors(analysis_callback2).l2_error[1, :], label = "with relaxation")
savefig("error_growth_relaxation.png") # hide
nothing # hide
```


![error growth relaxation](error_growth_relaxation.png)

For additional information, see [Ranocha et al. (2020)](https://doi.org/10.1137/19M1263480).

## References

[^Ketcheson2019]:
    Ketcheson (2019):
    Relaxation Runge-Kutta Methods: Conservation and stability for Inner-Product Norms.
    [DOI: 10.1137/19M1263662](https://doi.org/10.1137/19M1263662)

[^RanochaSayyariDalcinParsaniKetcheson2020]:
    Ranocha, Sayyari, Dalcin, Parsani, Ketcheson (2020):
    Relaxation Runge–Kutta Methods: Fully-Discrete Explicit Entropy-Stable Schemes for the Compressible Euler and Navier–Stokes Equations.
    [DOI: 10.1137/19M1263480](https://doi.org/10.1137/19M1263480)



# Plotting the solution

## conversion functions

when plotting can also provide a `conversion` function that converts the solution. A conversion function should take the values
of the primitive variables `q` at one node, and the `equations` as input and should return an `SVector` of any length as output. For a user defined conversion function,
there should also exist a function `varnames(conversion, equations)` that returns a `Tuple` of the variable names used for labelling. The conversion function can, e.g.,
be [`prim2cons`](@ref) or [`waterheight_total`](@ref) if one only wants to plot the total water height. The resulting plot will have one subplot for each of the returned
variables of the conversion variable. By default, the conversion function is just [`prim2phys`](@ref), which computes the physical variables
from the primitive ones. For most equations this is the identity, but for hyperbolic approximations like [`HyperbolicSerreGreenNaghdiEquations1D`](@ref), it returns a reduced set of variables (excluding auxiliary variables like `w` and `H`).
