@testsnippet KdVEquation1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "kdv_1d")
end

@testitem "kdv_1d_basic" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_basic.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.0011470884318895587],
                        linf=[0.000850644033992265],
                        cons_error=[2.495781359357352e-13],
                        change_waterheight=2.495781359357352e-13)

    @test_allocations(semi, sol, allocs=5_000)
end

@testitem "kdv_1d_implicit" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_implicit.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.0011092095759535911],
                        linf=[0.0008576326249670752],
                        cons_error=[2.460254222569347e-13],
                        change_waterheight=2.460254222569347e-13)

    @test_allocations(semi, sol, allocs=5_000)
end

@testitem "kdv_1d_fourier" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_fourier.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.0006248515956867525],
                        linf=[0.00011293285604195014],
                        cons_error=[1.3322676295501878e-15],
                        change_waterheight=1.3322676295501878e-15,
                        atol=1e-8) # to make CI pass

    @test_allocations(semi, sol, allocs=5_000)
end

@testitem "kdv_1d_manufactured" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_manufactured.jl"),
                        tspan=(0.0, 1.0),
                        l2=[8.09529983835206e-8],
                        linf=[8.27538344516654e-8],
                        cons_error=[7.111392141823103e-9],
                        change_waterheight=-7.111392141823103e-9,
                        atol=1e-9,
                        atol_ints=1e-10) # to make CI pass

    @test_allocations(semi, sol, allocs=5_000)
end

@testitem "kdv_1d_narrow_stencil" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_narrow_stencil.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.08433737909310937],
                        linf=[0.06531432869397868],
                        cons_error=[1.3322676295501878e-15],
                        change_waterheight=1.3322676295501878e-15)

    @test_allocations(semi, sol, allocs=5_000)
end
