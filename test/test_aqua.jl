@testitem "Aqua.jl" setup=[Setup] begin
    using Aqua
    using ExplicitImports: test_explicit_imports

    # Aqua.jl
    Aqua.test_all(DispersiveShallowWater,
                  ambiguities = false)

    # ExplicitImports.jl
    # Non-public names from SummationByPartsOperators.jl, which we need to dispatch
    # on the operator types and to extend `xmin` and `xmax`
    non_public_explicit_imports = (:AbstractDerivativeOperator,
                                   :AbstractNonperiodicDerivativeOperator,
                                   :AbstractPeriodicDerivativeOperator,
                                   :UniformCoupledOperator, :UniformPeriodicCoupledOperator,
                                   :xmin, :xmax)
    non_public_qualified_accesses = (
                                     # Base, used in `@autoinfiltrate`
                                     :PkgId, :eval, :loaded_modules,
                                     # Base, used in the `AnalysisCallback`
                                     :gc_live_bytes, :gc_time_ns,
                                     # Base, extended for broadcasting
                                     :broadcastable,
                                     # TimerOutputs.jl
                                     :complement!,
                                     # TrixiBase.jl
                                     :find_assignment,
                                     # ForwardDiff.jl
                                     :jacobian, :pickchunksize)
    # Names can only be declared `public` since Julia v1.11. Thus, we only check
    # whether the imported and accessed names are public on Julia v1.11 and newer.
    check_public = VERSION >= v"1.11"
    imports_are_public = check_public && (ignore = non_public_explicit_imports,)
    accesses_are_public = check_public && (ignore = non_public_qualified_accesses,)
    test_explicit_imports(DispersiveShallowWater;
                          no_stale_explicit_imports = (ignore = (:derivative_discontinuity!,),),
                          all_explicit_imports_are_public = imports_are_public,
                          all_qualified_accesses_are_public = accesses_are_public)
end
