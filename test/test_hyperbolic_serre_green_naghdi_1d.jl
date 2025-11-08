@testsnippet HyperbolicSerreGreenNaghdiEquations1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "hyperbolic_serre_green_naghdi_1d")
end

@testitem "hyperbolic_serre_green_naghdi_soliton.jl" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_soliton.jl"),
                        tspan=(0.0, 0.1),
                        l2=[
                            0.0007038714283663042,
                            0.006508261273448058,
                            0.0,
                            0.024517136865274798,
                            0.002141907410685252
                        ],
                        linf=[
                            0.0005088046605401519,
                            0.0036954890877776703,
                            0.0,
                            0.015022422297545818,
                            0.0013290414555349184
                        ],
                        cons_error=[
                            2.7000623958883807e-13,
                            0.00013389587974454997,
                            0.0,
                            0.005963937086921899,
                            4.502801745331908e-5
                        ],
                        change_entropy_modified=-2.3374946067633573e-7)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_soliton.jl with bathymetry_mild_slope" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_soliton.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        tspan=(0.0, 0.1),
                        l2=[
                            0.0007038714283663042,
                            0.006508261273448058,
                            0.0,
                            0.024517136865274798,
                            0.002141907410685252
                        ],
                        linf=[
                            0.0005088046605401519,
                            0.0036954890877776703,
                            0.0,
                            0.015022422297545818,
                            0.0013290414555349184
                        ],
                        cons_error=[
                            2.7000623958883807e-13,
                            0.00013389587974454997,
                            0.0,
                            0.005963937086921899,
                            4.502801745331908e-5
                        ],
                        change_entropy_modified=-2.3374946067633573e-7)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_soliton_relaxation.jl" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_soliton_relaxation.jl"),
                        tspan=(0.0, 0.1),
                        l2=[
                            0.0007041797674417557,
                            0.006510737539763134,
                            0.0,
                            0.024517447804525746,
                            0.002141928791106223
                        ],
                        linf=[
                            0.0005090662088376163,
                            0.0036987746989370907,
                            0.0,
                            0.01502004552088677,
                            0.0013289272946777064
                        ],
                        cons_error=[
                            3.126388037344441e-13,
                            0.00013409338344283483,
                            0.0,
                            0.005963706457799891,
                            4.504121848469822e-5
                        ],
                        change_entropy_modified=-5.684341886080802e-14,
                        atol=2.0e-10) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_soliton.jl Jacobian" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D,
    AdditionalImports
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_soliton.jl"),
                        tspan=(0.0, 0.1),)

    J = @test_nowarn DispersiveShallowWater.jacobian(semi)
    @test size(J, 1) == 5 * nnodes(semi.mesh)
    @test size(J, 2) == 5 * nnodes(semi.mesh)

    # Check whether the Jacobian agrees with a finite difference approximation
    ode = semidiscretize(semi, (0.0, 0.1))
    q = ode.u0
    dq = similar(q)
    DispersiveShallowWater.rhs!(dq, q, semi, 0.0)
    h = similar(q)
    sqrt_eps = sqrt(eps())
    for i in eachindex(h)
        h[i] = sqrt_eps * (rand() - 0.5)
        q[i] += h[i]
    end
    dq_plus_h = similar(q)
    DispersiveShallowWater.rhs!(dq_plus_h, q, semi, 0.0)
    @test maximum(abs, vec(dq) + J * vec(h) - vec(dq_plus_h)) < 1.0e-11
end

@testitem "hyperbolic_serre_green_naghdi_soliton_reflecting.jl" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_soliton_reflecting.jl"),
                        l2=[
                            0.060131215257107934,
                            2.298340621897689,
                            0.0,
                            0.7453781251158317,
                            0.05994684118460208
                        ],
                        linf=[
                            0.02694706083045495,
                            1.1249101579283711,
                            0.0,
                            0.3234791689092935,
                            0.02684912262878769
                        ],
                        cons_error=[
                            5.381650680646999e-11,
                            6.870195197880939,
                            0.0,
                            0.002803316321279915,
                            9.906302523177146e-6
                        ],
                        change_waterheight=-5.381650680646999e-11,
                        change_entropy_modified=-7.304469784230605e-5)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_soliton_reflecting_relaxation.jl with bathymetry_mild_slope" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_soliton_reflecting_relaxation.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        l2=[
                            0.060123240377457116,
                            2.2983598045973,
                            0.0,
                            0.7454649703743562,
                            0.05993756877756754
                        ],
                        linf=[
                            0.02694677000403889,
                            1.1249124192485027,
                            0.0,
                            0.323562217380842,
                            0.026847036583736994
                        ],
                        cons_error=[
                            5.397282620833721e-11,
                            6.870200683471968,
                            0.0,
                            0.0028274466604412555,
                            9.11803174119541e-6
                        ],
                        change_waterheight=-5.397282620833721e-11,
                        change_entropy_modified=5.684341886080802e-14)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_well_balanced.jl" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_well_balanced.jl"),
                        tspan=(0.0, 0.01),
                        l2=[
                            1.493401506906233e-14,
                            4.5636580728083475e-14,
                            4.217312212200215e-15,
                            5.399527148467818e-14,
                            4.646551952637425e-15
                        ],
                        linf=[
                            3.175237850427948e-14,
                            1.2230280990924922e-13,
                            6.661338147750939e-15,
                            1.2201879375620406e-13,
                            2.4646951146678475e-14
                        ],
                        cons_error=[
                            1.509903313490213e-14,
                            1.7256186536178874e-16,
                            2.220446049250313e-15,
                            6.788390857160712e-14,
                            2.220446049250313e-15
                        ],
                        change_entropy_modified=-3.197442310920451e-14,
                        lake_at_rest=1.833689450281284e-14)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_manufactured.jl" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_manufactured.jl"),
                        tspan=(0.0, 0.1),
                        l2=[
                            0.0003106551093994142,
                            5.369678232967609e-5,
                            0.0,
                            0.005791457102000912,
                            0.00031799455343148314
                        ],
                        linf=[
                            0.0005248180367165567,
                            0.00011353070870012694,
                            0.0,
                            0.010051964701901284,
                            0.0005326971020860327
                        ],
                        cons_error=[4.440892098500626e-16,
                            8.350198113236118e-6,
                            0.0,
                            1.141036530464996,
                            4.819814592771365e-6
                        ], atol=1e-8) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_manufactured_reflecting.jl with bathymetry_mild_slope" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_manufactured_reflecting.jl"),
                        bathymetry_type=bathymetry_mild_slope,
                        l2=[
                            0.0004804670128194017,
                            9.780609099451367e-5,
                            0.0,
                            1.128004060442867,
                            0.0004874625171967133
                        ],
                        linf=[
                            0.007938551537243654,
                            0.00015994533252292054,
                            0.0,
                            1.7376848592488843,
                            0.008062275831221655
                        ],
                        cons_error=[
                            4.999999985935993,
                            0.31828162626441847,
                            0.0,
                            0.020761884510929786,
                            5.000009387010711
                        ],
                        atol=1e-9) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_manufactured_reflecting.jl with bathymetry_flat" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_manufactured_reflecting.jl"),
                        bathymetry_type=bathymetry_flat,
                        l2=[
                            0.000308924023550926,
                            3.5906811918117744e-5,
                            0.0,
                            0.001453438157422651,
                            0.00031356077105954764
                        ],
                        linf=[
                            0.004819147897776155,
                            7.19964755682076e-5,
                            0.0,
                            0.025612277450607124,
                            0.0048963040276222
                        ],
                        cons_error=[
                            4.9999999859359985,
                            0.31829262810216014,
                            0.0,
                            0.9341006360094346,
                            5.000006253428024
                        ],
                        atol=1e-9) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_dingemans.jl" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_dingemans.jl"),
                        tspan=(0.0, 1.0),
                        l2=[
                            0.26477322527867336,
                            0.8633746824236788,
                            5.2597559566212386e-15,
                            0.5966739911067777,
                            0.26574526372150303
                        ],
                        linf=[
                            0.03634214058318297,
                            0.1179926788229248,
                            5.662137425588298e-15,
                            0.0803729516282009,
                            0.03647842275648405
                        ],
                        cons_error=[
                            2.7569058147491887e-12,
                            0.000940711448287545,
                            0.0,
                            0.048526294267482256,
                            0.00024128430453629335
                        ],
                        change_entropy=-0.0018582545795879923,
                        change_entropy_modified=-3.388064897080767e-6)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_serre_green_naghdi_conservation.jl" setup=[
    Setup,
    HyperbolicSerreGreenNaghdiEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_serre_green_naghdi_conservation.jl"),
                        l2=[
                            1.3601936429761978,
                            2.3678119582029846,
                            3.565537895102537e-14,
                            0.8095117331716646,
                            1.3613672882028665
                        ],
                        linf=[
                            1.0010267421677255,
                            0.7870064257071532,
                            9.325873406851315e-15,
                            0.24155794880613263,
                            1.0010070893015126
                        ],
                        cons_error=[6.161826604511589e-11,
                            0.0003412183411564129,
                            2.2737367544323206e-13,
                            0.007621698018910728,
                            0.0014365872115718048],
                        change_entropy=-0.20915660108039447,
                        change_entropy_modified=-0.09432092643828582,
                        atol=1e-8) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end
