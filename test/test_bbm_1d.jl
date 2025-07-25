@testsnippet BBMEquation1D begin
    EXAMPLES_DIR = joinpath(examples_dir(), "bbm_1d")
end

@testitem "bbm_1d_basic" setup=[Setup, BBMEquation1D, AdditionalImports] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_basic.jl"),
                        tspan=(0.0, 100.0),
                        l2=[0.020905920116663004],
                        linf=[0.011361008348288737],
                        cons_error=[1.3322676295501878e-15],
                        change_waterheight=1.3322676295501878e-15,
                        change_entropy_modified=-8.852448858398532e-7,
                        change_hamiltonian=-4.366246803000351e-6)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)

    # test upwind operators
    D1 = upwind_operators(periodic_derivative_operator; derivative_order = 1,
                          accuracy_order = accuracy_order, xmin = mesh.xmin,
                          xmax = mesh.xmax,
                          N = mesh.N)
    D2 = sparse(D1.minus) * sparse(D1.plus)
    solver = Solver(D1, D2)
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_basic.jl"),
                        tspan=(0.0, 100.0),
                        solver=solver,
                        l2=[0.11730278305145693],
                        linf=[0.06433115916307008],
                        cons_error=[6.306066779870889e-14],
                        change_waterheight=-6.306066779870889e-14,
                        change_entropy_modified=-0.005460856586602614,
                        change_hamiltonian=-0.0204980283782078)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)

    # test PeriodicRationalDerivativeOperator
    D1 = periodic_derivative_operator(1, accuracy_order, xmin(mesh), xmax(mesh),
                                      nnodes(mesh))
    D2 = D1^2
    solver = Solver(D1, D2)
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_basic.jl"),
                        tspan=(0.0, 100.0),
                        solver=solver,
                        l2=[0.013447144236536044],
                        linf=[0.007184057857459125],
                        cons_error=[1.5543122344752192e-15],
                        change_waterheight=-1.5543122344752192e-15,
                        change_entropy_modified=-8.68871272874383e-7,
                        change_hamiltonian=-3.697687313897191e-6)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)
end

@testitem "bbm_1d_basic with split_form = false" setup=[
    Setup,
    BBMEquation1D,
    AdditionalImports
] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_basic.jl"),
                        tspan=(0.0, 100.0),
                        split_form=false,
                        l2=[0.018773238618175637],
                        linf=[0.0102691348638429],
                        cons_error=[4.218847493575595e-15],
                        change_waterheight=4.218847493575595e-15,
                        change_entropy_modified=-5.764463673341158e-7,
                        change_hamiltonian=-3.3263833076890847e-6)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)

    # test upwind operators
    D1 = upwind_operators(periodic_derivative_operator; derivative_order = 1,
                          accuracy_order = accuracy_order, xmin = mesh.xmin,
                          xmax = mesh.xmax,
                          N = mesh.N)
    D2 = sparse(D1.minus) * sparse(D1.plus)
    solver = Solver(D1, D2)
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_basic.jl"),
                        tspan=(0.0, 100.0),
                        solver=solver,
                        split_form=false,
                        l2=[0.11994249550267427],
                        linf=[0.06576332657307044],
                        cons_error=[5.773159728050814e-14],
                        change_waterheight=-5.773159728050814e-14,
                        change_entropy_modified=-0.00545433375143814,
                        change_hamiltonian=-0.020473704519129887)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)
end

@testitem "bbm_1d_relaxation" setup=[Setup, BBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_relaxation.jl"),
                        tspan=(0.0, 100.0),
                        l2=[0.020884321819757214],
                        linf=[0.011350059432479132],
                        cons_error=[2.220446049250313e-15],
                        change_waterheight=2.220446049250313e-15,
                        change_entropy_modified=2.7755575615628914e-17,
                        change_hamiltonian=-1.039007974656947e-6)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)
end

@testitem "bbm_1d_hamiltonian_relaxation" setup=[Setup, BBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_hamiltonian_relaxation.jl"),
                        tspan=(0.0, 100.0),
                        l2=[0.018751641621494324],
                        linf=[0.010258179709480952],
                        cons_error=[3.552713678800501e-15],
                        change_waterheight=3.552713678800501e-15,
                        change_entropy_modified=3.085709884831367e-7,
                        change_hamiltonian=4.440892098500626e-16)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)
end

@testitem "bbm_1d_fourier" setup=[Setup, BBMEquation1D] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_fourier.jl"),
                        tspan=(0.0, 100.0),
                        l2=[1.812306511686721e-5],
                        linf=[1.0044426771715909e-5],
                        cons_error=[6.439293542825908e-15],
                        change_waterheight=-6.439293542825908e-15,
                        change_entropy_modified=-8.734545181576792e-7,
                        change_hamiltonian=-3.282887288680314e-6)

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_000)
end

@testitem "bbm_1d_manufactured" setup=[Setup, BBMEquation1D, AdditionalImports] begin
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_manufactured.jl"),
                        tspan=(0.0, 1.0),
                        l2=[6.073756646372025e-9],
                        linf=[8.869951328982495e-9],
                        cons_error=[2.60491292165127e-12],
                        change_waterheight=2.60491292165127e-12,
                        atol=1e-11) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=5_100)

    # test upwind operators
    D1 = upwind_operators(periodic_derivative_operator; derivative_order = 1,
                          accuracy_order = accuracy_order, xmin = mesh.xmin,
                          xmax = mesh.xmax,
                          N = mesh.N)
    D2 = sparse(D1.minus) * sparse(D1.plus)
    solver = Solver(D1, D2)
    @test_trixi_include(joinpath(EXAMPLES_DIR, "bbm_1d_manufactured.jl"),
                        tspan=(0.0, 1.0),
                        solver=solver,
                        l2=[9.94432411102018e-8],
                        linf=[1.026056097863659e-7],
                        cons_error=[3.745851908818973e-12],
                        change_waterheight=-3.745851908818973e-12,
                        atol=1e-11) # to make CI pass

    @test_allocations(DispersiveShallowWater.rhs!, semi, sol, allocs=6_000)
end
