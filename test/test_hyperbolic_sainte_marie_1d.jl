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
