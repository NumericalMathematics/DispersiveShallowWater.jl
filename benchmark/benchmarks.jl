using BenchmarkTools
using DispersiveShallowWater

const SUITE = BenchmarkGroup()

elixirs = [joinpath(examples_dir(), "bbm_1d", "bbm_1d_basic.jl"),
    joinpath(examples_dir(), "bbm_1d", "bbm_1d_fourier.jl"),
    joinpath(examples_dir(), "bbm_bbm_1d", "bbm_bbm_1d_dg.jl"),
    joinpath(examples_dir(), "bbm_bbm_1d", "bbm_bbm_1d_relaxation.jl"),
    joinpath(examples_dir(), "bbm_bbm_1d", "bbm_bbm_1d_upwind_relaxation.jl"),
    joinpath(examples_dir(), "bbm_bbm_1d", "bbm_bbm_1d_basic_reflecting.jl"),
    joinpath(examples_dir(), "hyperbolic_serre_green_naghdi_1d",
             "hyperbolic_serre_green_naghdi_dingemans.jl"),
    joinpath(examples_dir(), "serre_green_naghdi_1d",
             "serre_green_naghdi_well_balanced.jl"),
    joinpath(examples_dir(), "svaerd_kalisch_1d",
             "svaerd_kalisch_1d_dingemans_relaxation.jl"),
    joinpath(examples_dir(), "kdv_1d", "kdv_1d_basic.jl"),
    joinpath(examples_dir(), "kdv_1d", "kdv_1d_implicit.jl")]

for elixir in elixirs
    benchname = joinpath(basename(dirname(elixir)), basename(elixir)) * " - rhs!:"
    println("Running $benchname...")
    redirect_stdout(devnull) do
        trixi_include(@__MODULE__, elixir, tspan = (0.0, 1e-11))
        return nothing
    end
    # Since Julia v1.12, the global variables defined by `trixi_include` are newer
    # than the world age this loop runs in, so we need to access them via `invokelatest`.
    sol_elixir = @invokelatest (@__MODULE__).sol
    semi_elixir = @invokelatest (@__MODULE__).semi
    t = first(sol_elixir.prob.tspan)
    SUITE[benchname] = @benchmarkable DispersiveShallowWater.rhs!($(similar(sol_elixir.u[end])),
                                                                  $(copy(sol_elixir.u[end])),
                                                                  $(semi_elixir), $t)
end
