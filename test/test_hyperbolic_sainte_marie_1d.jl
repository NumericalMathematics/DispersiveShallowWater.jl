@testsnippet HyperbolicSainteMarieEquations1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "hyperbolic_sainte_marie_1d")
end

@testitem "hyperbolic_sainte_marie_conservation.jl" setup=[
    Setup,
    HyperbolicSainteMarieEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_sainte_marie_conservation.jl"),
                        tol=1.0e-12,
                        tspan=(0.0, 1.0),
                        l2=[
                            1.423280412274356,
                            2.3685775546329264,
                            7.26695561732344e-15,
                            0.7778533080445147,
                            3.4060674733738088
                        ],
                        linf=[
                            1.1050866054400976,
                            1.150984157070933,
                            2.1371793224034263e-15,
                            0.47878548301312557,
                            1.8555275982418526
                        ],
                        cons_error=[
                            8.560050446249079e-10,
                            2.8723044393164088e-5,
                            1.6986412276764894e-15,
                            0.4680239632418549,
                            8.59621928160862
                        ],
                        change_waterheight=-8.560050446249079e-10,
                        change_momentum=-0.000768001258112605,
                        change_entropy=-0.36374450470134434,
                        change_entropy_modified=-8.867345968610607e-9)

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
                            1.488406724541389,
                            2.024017579537045,
                            0.0,
                            0.8522342446876735,
                            4.0229842125125765
                        ],
                        linf=[
                            1.1623280008399592,
                            0.9969234037889885,
                            0.0,
                            0.3984623461616125,
                            1.799637261193249
                        ],
                        cons_error=[
                            8.105871529551223e-11,
                            2.0319512388233818e-5,
                            0.0,
                            1.1009759246326063,
                            7.8970857389160445
                        ],
                        change_waterheight=-8.105871529551223e-11,
                        change_momentum=-1.7186252421197423e-12,
                        change_entropy=-0.5378173222520672,
                        change_entropy_modified=-1.3169483281672e-9)

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
                            0.26643241429880193,
                            0.8758998593861568,
                            5.277989040277512e-15,
                            0.2977904521550398,
                            0.1324767069050882
                        ],
                        linf=[
                            0.03652784253085728,
                            0.11943351601871774,
                            5.662137425588298e-15,
                            0.04015026702892366,
                            0.03575062143105756
                        ],
                        cons_error=[
                            3.126388037344441e-12,
                            0.00028199582134299866,
                            0.0,
                            0.01651557336722678,
                            0.15898659132614137
                        ],
                        change_waterheight=-3.126388037344441e-12,
                        change_momentum=-3.527300661276822e-8,
                        change_entropy=-0.0007378805827329415,
                        change_entropy_modified=-7.357333515756181e-7)

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
                            0.00019537168346170865,
                            3.453172345767322e-5,
                            0.0,
                            0.002547202021354868,
                            0.001859407766133216
                        ],
                        linf=[
                            0.0003058377491629294,
                            6.150206698318783e-5,
                            0.0,
                            0.004072090053441002,
                            0.0030164198069354553
                        ],
                        cons_error=[
                            4.440892098500626e-16,
                            1.4747533280490277e-5,
                            0.0,
                            0.4936130523192781,
                            0.0001957021433037043
                        ])

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end

@testitem "hyperbolic_sainte_marie_manufactured.jl with bathymetry_flat" setup=[
    Setup,
    HyperbolicSainteMarieEquations1D
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR,
                                 "hyperbolic_sainte_marie_manufactured.jl"),
                        bathymetry_type=bathymetry_flat,
                        tspan=(0.0, 0.1),
                        abstol=1.0e-12,
                        reltol=1.0e-12,
                        l2=[
                            5.523682478151583e-5,
                            2.2904928037377655e-5,
                            0.0,
                            0.0007163282391786536,
                            0.0005017750313044287
                        ],
                        linf=[
                            9.87142981752065e-5,
                            3.7839638703074385e-5,
                            0.0,
                            0.0013307590180882123,
                            0.0009464535129504181
                        ],
                        cons_error=[
                            0.0,
                            2.567887582648476e-6,
                            0.0,
                            0.6474616401127097,
                            6.388777581955773e-5
                        ])

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=1_000)
end
