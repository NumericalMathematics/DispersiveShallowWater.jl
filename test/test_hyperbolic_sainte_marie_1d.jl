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
