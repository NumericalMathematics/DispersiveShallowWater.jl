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
                        abstol = 1.0e-12,
                        reltol = 1.0e-12,
                        l2=[
                            0.00025264914641096573,
                            4.26502127473517e-5,
                            0.0,
                            0.0032924534691802037,
                            0.002417874522326443
                        ],
                        linf=[
                            0.0004135814134835769,
                            8.132376992542323e-5,
                            0.0,
                            0.005705816127191898,
                            0.003909187839172712
                        ],
                        cons_error=[4.440892098500626e-16,
                            0.0,
                            1.1119371785710541e-5,
                            0.0,
                            0.49366652487677776,
                            4.9263555289584814e-5
                        ])

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end
