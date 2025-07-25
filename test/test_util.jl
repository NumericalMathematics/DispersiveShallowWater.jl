using Test: @test
using TrixiTest: @trixi_test_nowarn, @test_allocations, get_kwarg, append_to_kwargs

# Use a macro to avoid world age issues when defining new initial conditions etc.
# inside an example.
macro test_trixi_include_base(example, args...)
    local additional_ignore_content = get_kwarg(args, :additional_ignore_content, Any[])
    local l2 = get_kwarg(args, :l2, nothing)
    local linf = get_kwarg(args, :linf, nothing)
    local cons_error = get_kwarg(args, :cons_error, nothing)
    local change_waterheight = get_kwarg(args, :change_waterheight, nothing)
    local change_velocity = get_kwarg(args, :change_velocity, nothing)
    local change_momentum = get_kwarg(args, :change_momentum, nothing)
    local change_entropy = get_kwarg(args, :change_entropy, nothing)
    local change_entropy_modified = get_kwarg(args, :change_entropy_modified, nothing)
    local change_hamiltonian = get_kwarg(args, :change_hamiltonian, nothing)
    local lake_at_rest = get_kwarg(args, :lake_at_rest, nothing)
    local atol = get_kwarg(args, :atol, 1e-12)
    local rtol = get_kwarg(args, :rtol, sqrt(eps()))
    local atol_ints = get_kwarg(args, :atol_ints, 1e-11)
    local rtol_ints = get_kwarg(args, :rtol_ints, sqrt(eps()))

    local kwargs = Pair{Symbol, Any}[]
    for arg in args
        if (arg.head == :(=) &&
            !(arg.args[1] in (:additional_ignore_content,
                              :l2, :linf, :cons_error, :change_waterheight,
                              :change_velocity, :change_momentum, :change_entropy,
                              :change_entropy_modified, :change_hamiltonian,
                              :lake_at_rest,
                              :atol, :rtol, :atol_ints, :rtol_ints)))
            push!(kwargs, Pair(arg.args...))
        end
    end

    quote
        println("═"^100)
        println($example)

        # evaluate examples in the scope of the module they're called from
        @trixi_test_nowarn trixi_include(@__MODULE__, $example; $kwargs...) $additional_ignore_content

        # if present, compare l2, linf and conservation errors against reference values
        if !isnothing($l2) || !isnothing($linf) || !isnothing($cons_error)
            errs = errors(analysis_callback)

            if !isnothing($l2)
                l2_measured = errs.l2_error[:, end]
                @test length($l2) == length(l2_measured)
                for (l2_expected, l2_actual) in zip($l2, l2_measured)
                    @test isapprox(l2_expected, l2_actual, atol = $atol, rtol = $rtol)
                end
            end

            if !isnothing($linf)
                linf_measured = errs.linf_error[:, end]
                @test length($linf) == length(linf_measured)
                for (linf_expected, linf_actual) in zip($linf, linf_measured)
                    @test isapprox(linf_expected, linf_actual, atol = $atol, rtol = $rtol)
                end
            end

            if !isnothing($cons_error)
                cons_error_measured = errs.conservation_error[:, end]
                @test length($cons_error) == length(cons_error_measured)
                for (conservation_error_expected, conservation_error_actual) in zip($cons_error,
                                                                                    cons_error_measured)
                    @test isapprox(conservation_error_expected, conservation_error_actual,
                                   atol = $atol, rtol = $rtol)
                end
            end
        end

        if !isnothing($change_waterheight) || !isnothing($change_velocity) ||
           !isnothing($change_momentum) ||
           !isnothing($change_entropy) || !isnothing($change_entropy_modified) ||
           !isnothing($change_hamiltonian) || !isnothing($lake_at_rest)
            ints = integrals(analysis_callback)

            if !isnothing($change_waterheight)
                waterheight_change_measured = ints.waterheight_total[end] -
                                              ints.waterheight_total[1]
                @test isapprox($change_waterheight, waterheight_change_measured,
                               atol = $atol_ints, rtol = $rtol_ints)
            end

            if !isnothing($change_velocity)
                velocity_change_measured = ints.velocity[end] - ints.velocity[1]
                @test isapprox($change_velocity, velocity_change_measured,
                               atol = $atol_ints,
                               rtol = $rtol_ints)
            end

            if !isnothing($change_momentum)
                momentum_change_measured = ints.momentum[end] - ints.momentum[1]
                @test isapprox($change_momentum, momentum_change_measured,
                               atol = $atol_ints,
                               rtol = $rtol_ints)
            end

            if !isnothing($change_entropy)
                entropy_change_measured = ints.entropy[end] - ints.entropy[1]
                @test isapprox($change_entropy, entropy_change_measured, atol = $atol_ints,
                               rtol = $rtol_ints)
            end

            if !isnothing($change_entropy_modified)
                entropy_modified_change_measured = ints.entropy_modified[end] -
                                                   ints.entropy_modified[1]
                @test isapprox($change_entropy_modified, entropy_modified_change_measured,
                               atol = $atol_ints,
                               rtol = $rtol_ints)
            end

            if !isnothing($change_hamiltonian)
                change_hamiltonian_change_measured = ints.hamiltonian[end] -
                                                     ints.hamiltonian[1]
                @test isapprox($change_hamiltonian, change_hamiltonian_change_measured,
                               atol = $atol_ints,
                               rtol = $rtol_ints)
            end

            if !isnothing($lake_at_rest)
                lake_at_rest_measured = ints.lake_at_rest_error[end]
                @test isapprox($lake_at_rest, lake_at_rest_measured, atol = $atol_ints,
                               rtol = $rtol_ints)
            end
        end
        println("═"^100)
    end
end

"""
    @test_trixi_include(example; l2=nothing, linf=nothing, cons_error=nothing
                                change_waterheight=nothing,
                                change_velocity=nothing,
                                change_entropy=nothing,
                                change_entropy_modified=nothing,
                                change_hamiltonian=nothing,
                                lake_at_rest=nothing,
                                atol=1e-12, rtol=sqrt(eps()),
                                atol_ints=1e-11, rtol_ints=sqrt(eps()))

Test by calling `trixi_include(example; parameters...)`.
By default, only the absence of error output is checked.
If `l2`, `linf` or `cons_error` are specified, in addition the resulting L2/Linf/conservation
errors are compared approximately against these reference values, using `atol, rtol`
as absolute/relative tolerance.
If `change_waterheight`, `change_velocity`, `change_momemtum`, `change_entropy`, `change_entropy_modified`,
`change_hamiltonian`, or `lake_at_rest` are specified, in addition the resulting changes of the different errors are
compared approximately against these reference values, using `atol_ints`, `rtol_ints` as absolute/relative tolerance.
"""
macro test_trixi_include(expr, args...)
    local add_to_additional_ignore_content = [
        r"┌ Warning: The still-water surface needs to be 0 for the BBM-BBM equations\n└ @ DispersiveShallowWater .*\n"
    ]
    args = append_to_kwargs(args, :additional_ignore_content,
                            add_to_additional_ignore_content)
    quote
        @test_trixi_include_base($(esc(expr)), $(args...))
    end
end
