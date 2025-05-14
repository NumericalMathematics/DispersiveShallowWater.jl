@testsnippet KdVEquation1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "kdv_1d")
end

@testitem "kdv_1d_basic" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_basic.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.0011470884316580313],
                        linf=[0.0008506440340584343],
                        cons_error=[2.220446049250313e-15],
                        change_waterheight=2.220446049250313e-15,)

    @test_allocations(semi, sol, allocs=5_000)
end

@testitem "kdv_1d_basic using ForwardDiff" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_implicit.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.0011092095757747904],
                        linf=[0.0008576326250430144],
                        cons_error=[1.7763568394002505e-15],
                        change_waterheight=2.220446049250313e-15,)

    @test_allocations(semi, sol, allocs=5_000)
end

@testitem "kdv_1d_fourier" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_fourier.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.0006248516017854692],
                        linf=[0.00011293320517307625],
                        cons_error=[4.440892098500626e-15],
                        change_waterheight=4.440892098500626e-15,)

    @test_allocations(semi, sol, allocs=5_000)
end

err = errors(analysis_callback)
int = integrals(analysis_callback)

@testitem "kdv_1d_manufactured" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_manufactured.jl"),
                        tspan=(0.0, 1.0),
                        l2=[8.19260012441373e-8],
                        linf=[8.736123491281944e-8],
                        cons_error=[1.545513450551539e-8],
                        change_waterheight=-1.545513450551539e-8,
                        atol=1e-9) # to make CI pass)

    @test_allocations(semi, sol, allocs=5_000)
end

@testitem "kdv_1d_wide_stencil" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_wide_stencil.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.08433737909307822],
                        linf=[0.06531432869388931],
                        cons_error=[1.3322676295501878e-15],
                        change_waterheight=1.3322676295501878e-15,)

    @test_allocations(semi, sol, allocs=5_000)
end
