@testitem "Visualization" setup=[Setup] begin
    using Plots
    custom_integral(q, equations) = q[1]^2
    DispersiveShallowWater.pretty_form_utf(::typeof(custom_integral)) = "∫η²"
    trixi_include(@__MODULE__, default_example(), tspan = (0.0, 1.0),
                  extra_analysis_integrals = (waterheight_total, custom_integral))
    @test_nowarn plot(semi => sol)
    @test_nowarn plot!(semi => sol, plot_initial = true)
    @test_nowarn plot!(semi => sol, plot_analytical = true)
    @test_nowarn plot(semi, sol, conversion = prim2cons, plot_bathymetry = false)
    @test_nowarn plot(semi => sol, 0.0)
    @test_nowarn plot(semi, sol, 0.0, conversion = prim2cons)
    @test_nowarn plot(analysis_callback)
    @test_nowarn plot(analysis_callback, what = (:errors,))
    @test_nowarn plot(analysis_callback, what = (:integrals, :errors))

    # Variables after the bathymetry need to be plotted, too
    trixi_include(@__MODULE__,
                  joinpath(examples_dir(), "hyperbolic_serre_green_naghdi_1d",
                           "hyperbolic_serre_green_naghdi_soliton.jl"),
                  tspan = (0.0, 0.1))
    for conversion in (prim2prim, prim2cons)
        names = filter(!in(("D", "b")),
                       collect(DispersiveShallowWater.varnames(conversion, equations)))
        for p in (plot(semi => sol; conversion, plot_bathymetry = false),
                  plot(semi => sol, 0.0; conversion))
            @test [series[:label] for series in p.series_list] == names
            # each variable is plotted in its own subplot
            @test allunique(series[:subplot] for series in p.series_list)
        end
    end
end
