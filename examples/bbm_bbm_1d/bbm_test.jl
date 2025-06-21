using SummationByPartsOperators: AbstractNonperiodicDerivativeOperator, AbstractPeriodicDerivativeOperator
equations = BBMBBMEquations1D(bathymetry_type = bathymetry_variable, gravity = 9.81)

initial_condition = initial_condition_convergence_test
boundary_conditions = boundary_condition_reflecting
# boundary_conditions = boundary_condition_periodic

coordinates_min = -130.0
coordinates_max = 20.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N )

solver = Solver(mesh, 4)
isa(solver.D1, AbstractPeriodicDerivativeOperator)
semi = Semidiscretization(mesh, equations, initial_condition, solver, boundary_conditions = boundary_conditions)

tspan = (0.0, 25.0)
ode = semidiscretize(semi, tspan)
analysis_callback = AnalysisCallback(semi; interval = 10,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 velocity, entropy),
                                     io = devnull)
callbacks = CallbackSet(analysis_callback)

saveat = range(tspan..., length = 100)
sol = solve(ode, Tsit5(), abstol = 1e-7, reltol = 1e-7,
            save_everystep = false, callback = callbacks, saveat = saveat)


using Plots

plot(semi => sol,  ylims = (-0.8, 0.1),
)         
plot(semi => sol)

anim = @animate for step in 1:length(sol.u)
    plot(semi => sol, plot_initial = true, conversion = waterheight_total, step = step, ylims = (-0.8, 0.1))
end
gif(anim, fps = 25)