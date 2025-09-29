@testitem "Aqua.jl" setup=[Setup] begin
    using Aqua
    using ExplicitImports: check_no_implicit_imports, check_no_stale_explicit_imports

    # Aqua.jl
    Aqua.test_all(DispersiveShallowWater,
                  ambiguities = false)

    # ExplicitImports.jl
    @test isnothing(check_no_implicit_imports(DispersiveShallowWater))
    @test isnothing(check_no_stale_explicit_imports(DispersiveShallowWater))
end
