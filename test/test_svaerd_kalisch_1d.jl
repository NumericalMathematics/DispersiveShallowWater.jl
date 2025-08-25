
@testsnippet SvaerdKalischEquations1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "svaerd_kalisch_1d")
end

@testitem "svaerd_kalisch_1d_manufactured" setup=[Setup, SvaerdKalischEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_manufactured.jl"),
                        tspan=(0.0, 0.1),
                        l2=[3.3755102050606554e-6 2.320896961343125e-7 0.0],
                        linf=[4.908886917176503e-6 3.888671399332466e-7 0.0],
                        cons_error=[2.42861286636753e-16 1.9224170696150768e-7 0.0],
                        change_waterheight=-2.42861286636753e-16,
                        change_entropy=0.1868146724821993,
                        atol=1e-9) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=90_000)
end

@testitem "svaerd_kalisch_1d_basic_reflecting" setup=[
    Setup,
    SvaerdKalischEquations1D,
    AdditionalImports
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "svaerd_kalisch_1d_basic_reflecting.jl"),
                        tspan=(0.0, 1.0),
                        abstol=1e-12,
                        reltol=1e-12, # this example is relatively unstable with higher tolerances
                        l2=[5.402315494643131e-6 6.70558305594986e-8 0.0],
                        linf=[9.787585531206844e-5 1.460266869784954e-7 0.0],
                        cons_error=[1.0484871583691058e-9 0.5469460930247998 0.0],
                        change_waterheight=1.0484871583691058e-9,
                        change_entropy_modified=459.90372362340514,
                        atol=1e-11) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=650_000)

    # test upwind operators
    D1 = upwind_operators(Mattsson2017; derivative_order = 1,
                          accuracy_order = accuracy_order, xmin = mesh.xmin,
                          xmax = mesh.xmax,
                          N = mesh.N)
    solver = Solver(D1)
    @test_trixi_include(joinpath(EXAMPLES_DIR, "svaerd_kalisch_1d_basic_reflecting.jl"),
                        tspan=(0.0, 1.0),
                        solver=solver,
                        abstol=1e-12,
                        reltol=1e-12, # this example is relatively unstable with higher tolerances
                        l2=[5.862278175937948e-6 4.11195454078554e-9 0.0],
                        linf=[3.135228725170691e-5 8.797787950237668e-8 0.0],
                        cons_error=[1.700425028441774e-9 0.5469460935005555 0.0],
                        change_waterheight=-1.700425028441774e-9,
                        change_entropy_modified=459.9037221442321,
                        atol=1e-10) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=650_000)
end

@testitem "svaerd_kalisch_1d_dingemans" setup=[
    Setup,
    SvaerdKalischEquations1D,
    AdditionalImports
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_dingemans.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.26685829773477077, 0.8800952098489548, 0.0],
                        linf=[0.03665249779728463, 0.12048532901685113, 0.0],
                        cons_error=[8.526512829121202e-14, 6.496724002244977e-5, 0.0],
                        change_waterheight=-8.526512829121202e-14,
                        change_entropy=-0.0003347684281607144,
                        change_entropy_modified=-8.815959517960437e-9)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=350_000)

    # test PeriodicRationalDerivativeOperator
    D1 = periodic_derivative_operator(1, accuracy_order, xmin(mesh), xmax(mesh),
                                      nnodes(mesh))
    D2 = D1^2
    solver = Solver(D1, D2)
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_dingemans.jl"),
                        tspan=(0.0, 1.0),
                        solver=solver,
                        l2=[0.2668587324501861, 0.8800918552821204, 0.0],
                        linf=[0.03665281298157985, 0.12092636417680247, 0.0],
                        cons_error=[1.1368683772161603e-13, 6.498269475070908e-5, 0.0],
                        change_waterheight=-1.1368683772161603e-13,
                        change_entropy=-0.00033396742981040006,
                        change_entropy_modified=-9.819814295042306e-9)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=350_000)
end

@testitem "svaerd_kalisch_1d_dingemans_cg" setup=[Setup, SvaerdKalischEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_dingemans_cg.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.2669280328828958, 0.8803010025365813, 0.0],
                        linf=[0.03678993282064258, 0.12097067363770743, 0.0],
                        cons_error=[1.7053025658242404e-13, 6.420678045948553e-5, 0.0],
                        change_waterheight=-1.7053025658242404e-13,
                        change_entropy=-0.0003283147891579574,
                        change_entropy_modified=-1.8083028408000246e-9)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=750_000)
end

@testitem "svaerd_kalisch_1d_dingemans_fourier" setup=[Setup, SvaerdKalischEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_dingemans_fourier.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.2668948282001354, 0.8801953155243658, 0.0],
                        linf=[0.03665936863729946, 0.12083125907877076, 0.0],
                        cons_error=[1.1368683772161603e-13, 6.61417489483421e-5, 0.0],
                        change_waterheight=-1.1368683772161603e-13,
                        change_entropy=-0.00033368459492066904,
                        change_entropy_modified=-7.291419024113566e-9)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=13_000_000)
end

@testitem "svaerd_kalisch_1d_dingemans_upwind" setup=[Setup, SvaerdKalischEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_dingemans_upwind.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.26694765330284975, 0.8803371688530466, 0.0],
                        linf=[0.036669553278771194, 0.1212379289727636, 0.0],
                        cons_error=[5.684341886080802e-14, 6.444827576975179e-5, 0.0],
                        change_waterheight=-5.684341886080802e-14,
                        change_entropy=-0.00032897021856115316,
                        change_entropy_modified=-9.04344688024139e-9)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=350_000)
end

@testitem "svaerd_kalisch_1d_dingemans_relaxation" setup=[Setup, SvaerdKalischEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_dingemans_relaxation.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.2668583011029348, 0.8800952209358568, 0.0],
                        linf=[0.03665249830112549, 0.12048533060761354, 0.0],
                        cons_error=[8.526512829121202e-14, 6.4967241631659e-5, 0.0],
                        change_waterheight=-8.526512829121202e-14,
                        change_entropy=-0.000334760128112066,
                        change_entropy_modified=0.0)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=350_000)
end

@testitem "svaerd_kalisch_1d_well_balanced" setup=[Setup, SvaerdKalischEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "svaerd_kalisch_1d_well_balanced.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.0 1.135448143093612e-14 0.0],
                        linf=[0.0 8.133477278069499e-15 0.0],
                        cons_error=[0.0 1.6056589579882354e-14 0.0],
                        change_waterheight=0.0,
                        change_momentum=1.5679986322667355e-14,
                        change_entropy=0.0,
                        lake_at_rest=0.0)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=150_000)
end
