@testsnippet HyperbolicSainteMarieEquations1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "hyperbolic_sainte_marie_1d")
end

@testitem "hyperbolic_sainte_marie_conservation.jl" setup=[
    Setup,
    HyperbolicSainteMarieEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_sainte_marie_conservation.jl"),
                        l2=[
                            1.3572594173842452,
                            2.3633832882065207,
                            3.565537895102537e-14,
                            0.4622560613787918,
                            0.6144540736230862
                        ],
                        linf=[
                            1.001198219224467,
                            0.877550505633545,
                            9.325873406851315e-15,
                            0.16699894605274684,
                            0.27183886706644345
                        ],
                        cons_error=[
                            5.0874859880423173e-11,
                            0.00044041183843157583,
                            4.547473508864641e-13,
                            0.021403100433988752,
                            0.028407093718690702
                        ],
                        change_entropy=-0.2400751754171324,
                        change_entropy_modified=-0.12689501571117034,
                        atol=1e-8) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_sainte_marie_conservation.jl with bathymetry_flat" setup=[
    Setup,
    HyperbolicSainteMarieEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_sainte_marie_conservation.jl"),
                        b=0.0,
                        bathymetry_type=bathymetry_flat,
                        tol=1.0e-12,
                        tspan=(0.0, 1.0),
                        l2=[
                            1.4883751342392053,
                            2.0240715489993337,
                            5.384295521873979e-14,
                            0.851823022172266,
                            4.021158410315962
                        ],
                        linf=[
                            1.1622352722043399,
                            0.996970504640431,
                            3.1086244689504383e-15,
                            0.3985082548812433,
                            1.798277246494875
                        ],
                        cons_error=[
                            8.105871529551223e-11,
                            2.036050215359353e-5,
                            1.7053025658242404e-13,
                            1.1013360423538034,
                            7.893957675920978
                        ],
                        change_waterheight_total=-8.105871529551223e-11,
                        change_momentum=-1.730615650785694e-12,
                        change_entropy=-0.5373666791674623,
                        change_entropy_modified=-1.317857822868973e-9)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_sainte_marie_dingemans.jl" setup=[
    Setup,
    HyperbolicSainteMarieEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_sainte_marie_dingemans.jl"),
                        tspan=(0.0, 1.0),
                        l2=[
                            0.2664324191809125,
                            0.8758998976555538,
                            5.277989040277512e-15,
                            0.29779046756774213,
                            0.13247677055242074
                        ],
                        linf=[
                            0.036527826201786406,
                            0.1194334643284359,
                            5.662137425588298e-15,
                            0.0401503049607073,
                            0.03575050307041049
                        ],
                        cons_error=[
                            3.126388037344441e-12,
                            0.00028199575312887545,
                            0.0,
                            0.016515613345694284,
                            0.15898674767011647
                        ],
                        change_entropy=-0.0007378821217116638,
                        change_entropy_modified=-7.357332378887804e-7)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_sainte_marie_manufactured.jl" setup=[
    Setup,
    HyperbolicSainteMarieEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_sainte_marie_manufactured.jl"),
                        tspan=(0.0, 0.1),
                        abstol=1.0e-12,
                        reltol=1.0e-12,
                        l2=[
                            0.0002526491112439893,
                            4.265020844341071e-5,
                            0.0,
                            0.003292453459222419,
                            0.00241787560897696
                        ],
                        linf=[
                            0.00041358127213175777,
                            8.13234414983599e-5,
                            0.0,
                            0.005705816045885159,
                            0.003909190784468386
                        ],
                        cons_error=[
                            0.0,
                            1.1119366118657048e-5,
                            0.0,
                            0.49366652485721385,
                            4.926381115504918e-5
                        ])

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end
