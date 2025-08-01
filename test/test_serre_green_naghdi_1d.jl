@testsnippet SerreGreenNaghdiEquations1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "serre_green_naghdi_1d")
end

@testitem "serre_green_naghdi_soliton.jl" setup=[Setup, SerreGreenNaghdiEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton.jl"),
                        tspan=(0.0, 0.1),
                        l2=[9.994998669268741e-7, 1.4703973445698635e-6, 0.0],
                        linf=[6.5496216650196e-7, 1.027617322124641e-6, 0.0],
                        cons_error=[0.0, 8.174581012099225e-10, 0.0],
                        change_waterheight=0.0,
                        change_entropy_modified=-3.1093350116861984e-11)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=650_000)
end

@testitem "serre_green_naghdi_soliton.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same values as serre_green_naghdi_soliton.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        tspan=(0.0, 0.1),
                        l2=[9.994998669268741e-7, 1.4703973445698635e-6, 0.0],
                        linf=[6.5496216650196e-7, 1.027617322124641e-6, 0.0],
                        cons_error=[0.0, 8.174581012099225e-10, 0.0],
                        change_waterheight=0.0,
                        change_entropy_modified=-3.1093350116861984e-11)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=800_000)
end

@testitem "serre_green_naghdi_soliton.jl with bathymetry_variable" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same values as serre_green_naghdi_soliton.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton.jl"),
                        bathymetry_type=bathymetry_variable,
                        tspan=(0.0, 0.1),
                        l2=[9.994998669268741e-7, 1.4703973445698635e-6, 0.0],
                        linf=[6.5496216650196e-7, 1.027617322124641e-6, 0.0],
                        cons_error=[0.0, 8.174581012099225e-10, 0.0],
                        change_waterheight=0.0,
                        change_entropy_modified=-3.1093350116861984e-11)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=800_000)
end

@testitem "serre_green_naghdi_soliton_fourier.jl" setup=[Setup, SerreGreenNaghdiEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_fourier.jl"),
                        tspan=(0.0, 0.1),
                        l2=[8.252225014546995e-8, 6.724994492548714e-7, 0.0],
                        linf=[2.672093302180656e-8, 9.642725156897014e-8, 0.0],
                        cons_error=[2.842170943040401e-14, 4.627409566637652e-13, 0.0],
                        change_waterheight=2.842170943040401e-14,
                        change_entropy_modified=-3.097966327914037e-11)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=450_000)
end

@testitem "serre_green_naghdi_soliton_fourier.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same values as serre_green_naghdi_soliton_fourier.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_fourier.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        tspan=(0.0, 0.1),
                        l2=[8.252225014546995e-8, 6.724994492548714e-7, 0.0],
                        linf=[2.672093302180656e-8, 9.642725156897014e-8, 0.0],
                        cons_error=[2.842170943040401e-14, 4.627409566637652e-13, 0.0],
                        change_waterheight=2.842170943040401e-14,
                        change_entropy_modified=-3.097966327914037e-11)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=850_000)
end

@testitem "serre_green_naghdi_soliton_fourier.jl with bathymetry_variable" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same values as serre_green_naghdi_soliton_fourier.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_fourier.jl"),
                        bathymetry_type=bathymetry_variable,
                        tspan=(0.0, 0.1),
                        l2=[8.252225014546995e-8, 6.724994492548714e-7, 0.0],
                        linf=[2.672093302180656e-8, 9.642725156897014e-8, 0.0],
                        cons_error=[2.842170943040401e-14, 4.627409566637652e-13, 0.0],
                        change_waterheight=2.842170943040401e-14,
                        change_entropy_modified=-3.097966327914037e-11)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=850_000)
end

@testitem "serre_green_naghdi_soliton_upwind.jl" setup=[Setup, SerreGreenNaghdiEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_upwind.jl"),
                        tspan=(0.0, 0.1),
                        l2=[1.4876412924488654e-6, 5.9888097605442856e-6, 0.0],
                        linf=[1.0863034516361836e-6, 4.105927902009476e-6, 0.0],
                        cons_error=[4.263256414560601e-14, 4.483030568991353e-8, 0.0],
                        change_waterheight=4.263256414560601e-14,
                        change_entropy_modified=-3.1036506698001176e-11)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=600_000)
end

@testitem "serre_green_naghdi_soliton_upwind.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same as serre_green_naghdi_soliton_upwind.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_upwind.jl"),
                        tspan=(0.0, 0.1),
                        bathymetry_type=bathymetry_mild_slope,
                        l2=[1.4876412924488654e-6, 5.9888097605442856e-6, 0.0],
                        linf=[1.0863034516361836e-6, 4.105927902009476e-6, 0.0],
                        cons_error=[4.263256414560601e-14, 4.483030568991353e-8, 0.0],
                        change_waterheight=4.263256414560601e-14,
                        change_entropy_modified=-3.1036506698001176e-11)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=750_000)
end

@testitem "serre_green_naghdi_soliton_upwind.jl with bathymetry_variable" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same as serre_green_naghdi_soliton_upwind.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_upwind.jl"),
                        tspan=(0.0, 0.1),
                        bathymetry_type=bathymetry_variable,
                        l2=[1.4876412924488654e-6, 5.9888097605442856e-6, 0.0],
                        linf=[1.0863034516361836e-6, 4.105927902009476e-6, 0.0],
                        cons_error=[4.263256414560601e-14, 4.483030568991353e-8, 0.0],
                        change_waterheight=4.263256414560601e-14,
                        change_entropy_modified=-3.1036506698001176e-11)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=750_000)
end

@testitem "serre_green_naghdi_soliton_relaxation.jl" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_relaxation.jl"),
                        tspan=(0.0, 0.1),
                        l2=[8.252225169608892e-8, 6.724994488577288e-7, 0.0],
                        linf=[2.6716495016287922e-8, 9.642466235713909e-8, 0.0],
                        cons_error=[2.842170943040401e-14, 4.649614027130156e-13, 0.0],
                        change_waterheight=2.842170943040401e-14,
                        change_entropy_modified=0.0)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=450_000)
end

@testitem "serre_green_naghdi_soliton_relaxation.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same as serre_green_naghdi_soliton_relaxation.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_relaxation.jl"),
                        tspan=(0.0, 0.1),
                        bathymetry_type=bathymetry_mild_slope,
                        l2=[8.252225169608892e-8, 6.724994488577288e-7, 0.0],
                        linf=[2.6716495016287922e-8, 9.642466235713909e-8, 0.0],
                        cons_error=[2.842170943040401e-14, 4.649614027130156e-13, 0.0],
                        change_waterheight=2.842170943040401e-14,
                        change_entropy_modified=0.0)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=850_000)
end

@testitem "serre_green_naghdi_soliton_relaxation.jl with bathymetry_variable" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    # same as serre_green_naghdi_soliton_relaxation.jl but with more allocations
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_relaxation.jl"),
                        tspan=(0.0, 0.1),
                        bathymetry_type=bathymetry_variable,
                        l2=[8.252225169608892e-8, 6.724994488577288e-7, 0.0],
                        linf=[2.6716495016287922e-8, 9.642466235713909e-8, 0.0],
                        cons_error=[2.842170943040401e-14, 4.649614027130156e-13, 0.0],
                        change_waterheight=2.842170943040401e-14,
                        change_entropy_modified=0.0)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=850_000)
end

@testitem "serre_green_naghdi_well_balanced.jl" setup=[Setup, SerreGreenNaghdiEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_well_balanced.jl"),
                        tspan=(0.0, 2.0),
                        l2=[0, 0, 0],
                        linf=[0, 0, 0],
                        cons_error=[0, 0, 0],
                        change_waterheight=0.0,
                        change_momentum=0.0,
                        change_entropy_modified=0.0,
                        lake_at_rest=0.0)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=750_000)
end

@testitem "serre_green_naghdi_well_balanced.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_well_balanced.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        tspan=(0.0, 2.0),
                        l2=[0, 0, 0],
                        linf=[0, 0, 0],
                        cons_error=[0, 0, 0],
                        change_waterheight=0.0,
                        change_momentum=0.0,
                        change_entropy_modified=0.0,
                        lake_at_rest=0.0)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=750_000)
end

@testitem "serre_green_naghdi_dingemans.jl" setup=[Setup, SerreGreenNaghdiEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_dingemans.jl"),
                        tspan=(0.0, 1.0),
                        l2=[0.22632930215585131, 0.7400070292134782, 0.0],
                        linf=[0.036351214376643126, 0.11899056101300992, 0.0],
                        cons_error=[1.4210854715202004e-13, 3.194346928167053e-5, 0.0],
                        change_waterheight=-1.4210854715202004e-13,
                        change_entropy=2.282635693973134e-5,
                        change_entropy_modified=-9.135646905633621e-9)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=750_000)
end

@testitem "serre_green_naghdi_dingemans.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_dingemans.jl"),
                        tspan=(0.0, 1.0),
                        bathymetry_type=bathymetry_mild_slope,
                        l2=[0.22632930215585131, 0.7400070292134782, 0.0],
                        linf=[0.036351214376643126, 0.11899056101300992, 0.0],
                        cons_error=[1.4210854715202004e-13, 3.194346928167053e-5, 0.0],
                        change_waterheight=-1.4210854715202004e-13,
                        change_entropy=2.282635693973134e-5,
                        change_entropy_modified=-9.135646905633621e-9)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=750_000)
end

@testitem "serre_green_naghdi_conservation.jl" setup=[Setup, SerreGreenNaghdiEquations1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_conservation.jl"),
                        l2=[1.3655498085989206, 2.3967486930606716, 0.0],
                        linf=[1.001076318001934, 0.8052527556023067, 0.0],
                        cons_error=[0.0, 0.0002674927404067162, 0.0],
                        change_entropy=-0.05841897226287074,
                        change_entropy_modified=0.059273551933074486,
                        atol_ints=2e-8, # to make CI pass
                        atol=2e-8) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=900_000)
end

@testitem "serre_green_naghdi_conservation.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_conservation.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        l2=[1.3655493671985637, 2.3967828251339003, 0.0],
                        linf=[1.001075913983051, 0.8052680970114169, 0.0],
                        cons_error=[1.1368683772161603e-13, 0.00026407261543415217, 0.0],
                        change_entropy=-0.058352273553509804,
                        change_entropy_modified=0.05927340849780194,
                        atol_ints=2e-8) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=900_000)
end

@testitem "serre_green_naghdi_manufactured.jl with bathymetry_flat" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "serre_green_naghdi_manufactured.jl"),
                        bathymetry_type=bathymetry_flat,
                        l2=[9.264487112500761e-7 3.2951669640300735e-7 0.0],
                        linf=[2.3378983966537703e-6 5.457361316185683e-7 0.0],
                        cons_error=[2.6645352591003757e-15 2.9640116527840377e-7 0.0],
                        change_waterheight=-2.6645352591003757e-15,
                        change_entropy_modified=-22.793274391960267,
                        atol=1e-9) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=300_000)
end

@testitem "serre_green_naghdi_manufactured.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "serre_green_naghdi_manufactured.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        l2=[0.0029791064498190267 0.0004564521150416161 0.0],
                        linf=[0.007014500883379959 0.000676248180478678 0.0],
                        cons_error=[1.7763568394002505e-15 0.00019895120429092028 0.0],
                        change_waterheight=-1.7763568394002505e-15,
                        change_entropy_modified=135.16210732695845,
                        atol=1e-9) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=410_000)
end

@testitem "serre_green_naghdi_manufactured.jl with bathymetry_variable" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "serre_green_naghdi_manufactured.jl"),
                        bathymetry_type=bathymetry_variable,
                        l2=[0.00020091567099272015 6.26512426049926e-5 0.0],
                        linf=[0.0004699841860116294 0.00011399758749819355 0.0],
                        cons_error=[1.7763568394002505e-15 5.265416110067367e-5 0.0],
                        change_waterheight=-1.7763568394002505e-15,
                        change_entropy_modified=130.79560136094597,
                        atol=1e-11) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=410_000)
end

@testitem "serre_green_naghdi_manufactured_reflecting.jl with bathymetry_flat" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_manufactured_reflecting.jl"),
                        bathymetry_type=bathymetry_flat,
                        l2=[0.027195981356871006, 0.0003750952088972156, 0.0],
                        linf=[0.2907059851398657, 0.0007619654145489541, 0.0],
                        change_waterheight=15.970779079947338,
                        change_entropy_modified=2415.6322361337952)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=35_000)
end

@testitem "serre_green_naghdi_manufactured_reflecting.jl with bathymetry_variable" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_manufactured_reflecting.jl"),
                        bathymetry_type=bathymetry_variable,
                        l2=[0.008682031284538465, 0.0016019788677570093, 0.0],
                        linf=[0.07998918729370352, 0.0034380472896646253, 0.0],
                        change_waterheight=5.000659665129857,
                        change_entropy_modified=275.22930126732166)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=120_000)
end

@testitem "serre_green_naghdi_manufactured_reflecting_upwind.jl with bathymetry_flat" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_manufactured_reflecting_upwind.jl"),
                        bathymetry_type=bathymetry_flat,
                        l2=[0.039321533516191214, 0.0003059380760440109, 0.0],
                        linf=[0.3259202369703047, 0.0005643244796372793, 0.0],
                        change_waterheight=15.976311354158419,
                        change_entropy_modified=2416.856672789222)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=50_000)
end

@testitem "serre_green_naghdi_manufactured_reflecting_upwind.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_manufactured_reflecting.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        l2=[0.008756092883075014, 0.0016195432021351244, 0.0],
                        linf=[0.0804388379348362, 0.0034558532136001463, 0.0],
                        change_waterheight=5.000659665129858,
                        change_entropy_modified=274.83428431803526)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=120_000)
end

@testitem "serre_green_naghdi_manufactured_reflecting_var_coef_op.jl with bathymetry_flat" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_manufactured_reflecting_var_coef_op.jl"),
                        bathymetry_type=bathymetry_flat,
                        l2=[0.047962081783820484, 0.0002741366926583432, 0.0],
                        linf=[0.37478770827387464, 0.0004746167686701802, 0.0])

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=10_000)
end

@testitem "serre_green_naghdi_soliton_reflecting.jl with bathymetry_flat" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_reflecting.jl"),
                        bathymetry_type=bathymetry_flat,
                        l2=[0.10292128993725594, 0.2889327496936634, 0.0],
                        linf=[0.04452445088946688, 0.12376294777800739, 0.0],
                        change_waterheight=1.4210854715202004e-14,
                        change_entropy_modified=-1.7991674781114853e-6)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=260_000)
end

@testitem "serre_green_naghdi_soliton_reflecting.jl with bathymetry_mild_slope" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_reflecting.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        l2=[0.10292128993725498, 0.2889327496936603, 0.0],
                        linf=[0.044524450889467326, 0.12376294777800284, 0.0],
                        change_waterheight=1.4210854715202004e-14,
                        change_entropy_modified=-1.7991674781114853e-6)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=400_000)
end

@testitem "serre_green_naghdi_soliton_reflecting_var_coef_op.jl with bathymetry_flat" setup=[
    Setup,
    SerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "serre_green_naghdi_soliton_reflecting_var_coef_op.jl"),
                        bathymetry_type=bathymetry_flat,
                        l2=[0.10552684163688769, 0.29652604909005625, 0.0],
                        linf=[0.04598348365474836, 0.12697474293455474, 0.0],
                        change_waterheight=4.263256414560601e-14,
                        change_entropy_modified=-1.8172034401686687e-6,)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=60_000)
end
