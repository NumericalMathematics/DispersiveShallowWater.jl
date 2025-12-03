@testsnippet KdVEquation1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "kdv_1d")
end

@testitem "kdv_1d_basic" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_basic.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.0007835879713461127],
                        linf=[0.0005961613764722262],
                        cons_error=[4.440892098500626e-16],
                        change_waterheight=-4.440892098500626e-16,
                        change_entropy=3.664210623810504e-8)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)
end

@testitem "kdv_1d_relaxation" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_relaxation.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.000783588283119431],
                        linf=[0.0005961579568657394],
                        cons_error=[1.3322676295501878e-15],
                        change_waterheight=-1.3322676295501878e-15,
                        change_entropy=8.881784197001252e-16)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)
end

@testitem "kdv_1d_implicit" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_implicit.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.0007767194275956272],
                        linf=[0.0005970865295001904],
                        cons_error=[0.0],
                        change_waterheight=0.0,
                        change_entropy=-2.9601316864358296e-7)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)
end

@testitem "kdv_1d_fourier" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_fourier.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.0006248515956867525],
                        linf=[0.00011293285604195014],
                        cons_error=[1.3322676295501878e-15],
                        change_waterheight=1.3322676295501878e-15,
                        change_entropy=1.915115131367884e-6,
                        atol=1e-8) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)
end

@testitem "kdv_1d_manufactured" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_manufactured.jl"),
                        tspan=(0.0, 1.0),
                        l2=[5.366538289821288e-8],
                        linf=[5.472506603432237e-8],
                        cons_error=[8.881784197001252e-16],
                        change_waterheight=8.881784197001252e-16,
                        atol=1e-9,
                        atol_ints=1e-10) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)
end

@testitem "kdv_1d_narrow_stencil" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_narrow_stencil.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.0029648248065358984],
                        linf=[0.0020800831847420653],
                        cons_error=[8.881784197001252e-16],
                        change_waterheight=8.881784197001252e-16,
                        change_entropy=1.116256171940222e-7)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)
end

@testitem "kdv_1d_non_dim" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_non_dim.jl"),
                        tspan=(0.0, 50.0),
                        l2=[0.00031784195202372095],
                        linf=[8.131770758672274e-5],
                        cons_error=[2.842170943040401e-14],
                        change_waterheight=2.842170943040401e-14,
                        change_entropy=5.861043206323302e-9,
                        atol=1e-9) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)
end

@testitem "kdv_1d_IMEX" setup=[Setup, KdVEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "kdv_1d_IMEX.jl"),
                        tspan=(0.0, 5.0),
                        l2=[0.004952174509850488],
                        linf=[0.003962890861977875],
                        cons_error=[2.220446049250313e-15],
                        change_waterheight=-2.220446049250313e-15)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)
end
