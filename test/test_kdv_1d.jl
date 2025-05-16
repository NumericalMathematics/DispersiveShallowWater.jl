@testsnippet KdVEquation1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "kdv_1d")
end

@testitem "kdv_1d_basic" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_basic.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.0007835879714086921],
                        linf=[0.0005961613764309259],
                        cons_error=[1.056932319443149e-13],
                        change_waterheight=-1.056932319443149e-13)

    @test_allocations(semi, sol, allocs=5_000)
end

@testitem "kdv_1d_implicit" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_implicit.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.0007767194276962612],
                        linf=[0.0005970865294682159],
                        cons_error=[1.1057821325266559e-13],
                        change_waterheight=-1.1057821325266559e-13)

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
                        l2=[5.365647127534639e-8],
                        linf=[5.476293374329089e-8],
                        cons_error=[1.0208189848981419e-10],
                        change_waterheight=-1.0208189848981419e-10,
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
