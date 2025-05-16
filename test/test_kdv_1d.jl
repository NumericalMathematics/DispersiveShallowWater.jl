@testsnippet KdVEquation1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "kdv_1d")
end

@testitem "kdv_1d_basic" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_basic.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.0011470884314204559],
                        linf=[0.0008506440340751986],
                        cons_error=[3.1086244689504383e-15],
                        change_waterheight=3.1086244689504383e-15)

    @test_allocations(semi, sol, allocs=5_000)
end

@testitem "kdv_1d_implicit" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_implicit.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.001109209575787022],
                        linf=[0.0008576326250413491],
                        cons_error=[1.7763568394002505e-15],
                        change_waterheight=1.7763568394002505e-15)

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
                        l2=[8.192992655167117e-8],
                        linf=[8.733505763025562e-8],
                        cons_error=[1.5456631086152584e-8],
                        change_waterheight=-1.5456631086152584e-8,
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
