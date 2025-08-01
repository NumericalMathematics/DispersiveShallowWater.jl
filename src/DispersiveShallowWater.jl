"""
    DispersiveShallowWater

**DispersiveShallowWater.jl** is a Julia package that implements structure-preserving numerical methods for dispersive shallow water models.
It provides provably conservative, entropy-conserving, and well-balanced numerical schemes for some dispersive shallow water models.

The semidiscretizations are based on summation-by-parts (SBP) operators, which are implemented in
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl). To
obtain fully discrete schemes, the time integration methods from
[OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl)
are used to solve the resulting ordinary differential equations.
Fully discrete entropy-conservative methods can be obtained by using the relaxation method provided by DispersiveShallowWater.jl.

See also: [DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl)
"""
module DispersiveShallowWater

using BandedMatrices: BandedMatrix
using DiffEqBase: DiffEqBase, terminate!
using FastBroadcast: @..
using ForwardDiff: ForwardDiff
using Interpolations: Interpolations, linear_interpolation
using LinearAlgebra: mul!, ldiv!, I, Diagonal, Symmetric, diag, lu, cholesky, cholesky!,
                     issuccess
using PolynomialBases: PolynomialBases
using Printf: @printf, @sprintf
using PreallocationTools: DiffCache, get_tmp
using RecipesBase: RecipesBase, @recipe, @series
using RecursiveArrayTools: ArrayPartition
using Reexport: @reexport
using Roots: AlefeldPotraShi, find_zero

using SciMLBase: SciMLBase, DiscreteCallback, ODEProblem, ODESolution
import SciMLBase: u_modified!

@reexport using StaticArrays: SVector
using SimpleUnPack: @unpack
using SparseArrays: sparse, issparse
using SummationByPartsOperators: SummationByPartsOperators,
                                 AbstractDerivativeOperator,
                                 AbstractPeriodicDerivativeOperator,
                                 AbstractNonperiodicDerivativeOperator,
                                 PeriodicDerivativeOperator, PeriodicUpwindOperators,
                                 UniformPeriodicCoupledOperator,
                                 DerivativeOperator, UpwindOperators,
                                 UniformCoupledOperator,
                                 FourierDerivativeOperator,
                                 VarCoefDerivativeOperator,
                                 periodic_derivative_operator,
                                 derivative_order, integrate, mass_matrix,
                                 scale_by_mass_matrix!,
                                 scale_by_inverse_mass_matrix!,
                                 left_boundary_weight, right_boundary_weight
import SummationByPartsOperators: grid, xmin, xmax, semidiscretize
using TimerOutputs: TimerOutputs, print_timer, reset_timer!
@reexport using TrixiBase: trixi_include
using TrixiBase: TrixiBase, @trixi_timeit, timer

include("boundary_conditions.jl")
include("mesh.jl")
include("equations/equations.jl")
include("dispersion_relation.jl")
include("solver.jl")
include("semidiscretization.jl")
include("callbacks_step/callbacks_step.jl")
include("visualization.jl")
include("util.jl")

export examples_dir, get_examples, default_example, convergence_test

export AbstractShallowWaterEquations,
       KdVEquation1D,
       BBMEquation1D, BBMBBMEquations1D,
       SvärdKalischEquations1D, SvaerdKalischEquations1D,
       SerreGreenNaghdiEquations1D, HyperbolicSerreGreenNaghdiEquations1D

export LinearDispersionRelation, EulerEquations1D, wave_speed

export prim2prim, prim2cons, cons2prim, prim2phys,
       prim2nondim, nondim2prim,
       waterheight_total, waterheight,
       velocity, momentum, discharge,
       gravity,
       bathymetry, still_water_surface,
       energy_total, entropy, lake_at_rest_error,
       energy_total_modified, entropy_modified,
       hamiltonian

export Mesh1D, xmin, xmax, nnodes

export Solver

export Semidiscretization, semidiscretize, grid

export boundary_condition_periodic, boundary_condition_reflecting

export bathymetry_flat, bathymetry_mild_slope, bathymetry_variable

export initial_condition_convergence_test,
       initial_condition_soliton,
       initial_condition_manufactured, source_terms_manufactured,
       initial_condition_manufactured_reflecting, source_terms_manufactured_reflecting,
       initial_condition_dingemans,
       initial_condition_discontinuous_well_balancedness

export AnalysisCallback, RelaxationCallback, SummaryCallback
export tstops, errors, integrals

end
