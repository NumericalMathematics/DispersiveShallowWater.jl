# [Solvers](@id solvers)

This chapter covers different solvers and how to use them in [DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl).

To learn more about the analytical and mathematical background, go to the chapter about [Summation-by-Parts Operators](@ref sbp_operators).

## Introduction

[DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl) uses [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl/) as the foundation for all spatial discretizations. This package provides a comprehensive collection of summation-by-parts (SBP) operators that enable structure-preserving numerical methods with provable stability and conservation properties.

The [`Solver`](@ref) struct in [DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl) wraps the SBP operators needed for spatial discretization: typically a first-derivative operator `D1`, and optionally second- and third-derivative operators `D2` and `D3`, depending on the equation system being solved.

## Design Philosophy: Flexibility and Modularity

[DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl) follows a modular design that separates the choice of spatial discretization from the underlying equation systems. This design allows users to:

- **Compare discretization methods**: Test finite difference, discontinuous Galerkin, continuous Galerkin, and Fourier spectral approaches on identical problems
- **Investigate operator properties**: Analyze how different SBP operators preserve conservation laws and stability properties

This modularity is particularly relevant for research applications where the choice of discretization can significantly impact the preservation of physical invariants and long-time numerical stability. The unified interface enables systematic studies of numerical method performance without requiring separate implementations for each approach.

## Basic Summation-by-Parts Solver

As we have already seen in the [basic example](@ref basic_example), the easiest way to create a solver is to pass both a mesh and the desired accuracy order to the [`Solver`](@ref) function. This creates first-, second-, and third-derivative periodic summation-by-parts operators of the given accuracy order on the specified mesh:

```julia
accuracy_order = 4
solver = Solver(mesh, accuracy_order)
```

This approach creates periodic SBP operators and **only works with periodic boundary conditions**. If a solver containing periodic SBP operators is passed to `Semidiscretization` with non-periodic boundary conditions, it will throw an error.

!!! note "Accuracy Order Selection"
    Even accuracy orders are generally preferred for central difference operators as they provide symmetric stencils. Odd accuracy orders may result in reduced convergence rates due to asymmetric stencils.

## [Customizing Solvers](@id customize_solver)

For more advanced use cases, you can create custom solvers by explicitly constructing the required SBP operators. The general pattern is:

```julia
solver = Solver(D1, D2, D3)
```

where:

- `D1` is always required and must be an `AbstractDerivativeOperator`
- `D2` and `D3` are optional and can be either `AbstractDerivativeOperator`s, `AbstractMatrix`es, or `nothing`

## Reflecting Boundary Conditions

For non-periodic boundary conditions, you need to use non-periodic SBP operators. Here's how to create a solver for reflecting boundary conditions:

```julia
using SummationByPartsOperators: MattssonNordström2004, derivative_operator

# Create solver with SBP operators of accuracy order 4
accuracy_order = 4
D1 = derivative_operator(MattssonNordström2004(),
                         derivative_order = 1, accuracy_order = accuracy_order,
                         xmin = mesh.xmin, xmax = mesh.xmax, N = mesh.N)
solver = Solver(D1)
```

The `MattssonNordström2004()` operator family provides SBP operators that satisfy the SBP property with non-periodic boundary conditions, making them suitable for reflecting boundary conditions.
Other possible choices for sources can be found in the [documentation of SummationByPartsOperators.jl](https://ranocha.de/SummationByPartsOperators.jl/stable/).

## Upwind Operators

For equations that benefit from upwind discretizations (such as the Serre-Green-Naghdi equations), you can use [upwind SBP Operators](@ref upwind_sbp):

**For periodic boundary conditions:**

```julia
using SummationByPartsOperators: upwind_operators, periodic_derivative_operator

accuracy_order = 4
D1 = upwind_operators(periodic_derivative_operator;
                      derivative_order = 1,
                      accuracy_order = accuracy_order,
                      xmin = xmin(mesh), xmax = xmax(mesh),
                      N = nnodes(mesh))
solver = Solver(D1)
```

**For reflecting boundary conditions:**

```julia
using SummationByPartsOperators: Mattsson2017, upwind_operators

accuracy_order = 2
D1 = upwind_operators(Mattsson2017, derivative_order = 1, accuracy_order = accuracy_order,
                      xmin = xmin(mesh), xmax = xmax(mesh),
                      N = nnodes(mesh))
solver = Solver(D1)
```

Upwind operators provide additional dissipation and often achieve higher resolution for sharp features compared to central operators.

## KdV Equations with Third-Derivative Operators

The KdV equation requires first- and third-derivative operators, but no second-derivative operator. You can specify this by passing `nothing` for the second argument:

```julia
using SummationByPartsOperators: fourier_derivative_operator

# Create solver with Fourier SBP operators
D1 = fourier_derivative_operator(xmin(mesh), xmax(mesh), nnodes(mesh))
D3 = D1^3  # Third derivative via composition

solver = Solver(D1, nothing, D3)
```

## Using Sparse Matrices for Derivative Operators

While `D1` must always be an SBP operator, `D2` and `D3` can be regular (sparse) matrices. This can be useful for creating custom discretizations:

```julia
using SummationByPartsOperators: upwind_operators, periodic_derivative_operator
using SparseArrays: sparse

accuracy_order = 4
D1 = upwind_operators(periodic_derivative_operator; derivative_order = 1,
                      accuracy_order = accuracy_order, xmin = mesh.xmin, xmax = mesh.xmax,
                      N = mesh.N)
D2 = sparse(D1.plus) * sparse(D1.minus)  # Create D2 as sparse matrix
solver = Solver(D1, D2)
```

!!! warning "Performance consideration"
    While using sparse matrices for `D2` and `D3` provides flexibility, SBP operator-vector products are typically about an order of magnitude faster than sparse matrix-vector products. See the [SummationByPartsOperators.jl benchmarks](https://ranocha.de/SummationByPartsOperators.jl/stable/benchmarks/) for details.
    This does not matter if the operators provided as matrices are not required for operator/matrix-vector products but only for linear solves, e.g., for the BBM-BBM equations.

## Discontinuous and continuous Galerkin Methods

[DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl) supports [discontinuous Galerkin (DG)](@ref dg_sbp) and [continuous Galerkin (CG)](@ref cg_sbp) methods through coupled Lobatto-Legendre operators:

```julia
using SummationByPartsOperators: legendre_derivative_operator,
                                 UniformPeriodicMesh1D,
                                 couple_discontinuously
using SparseArrays: sparse

# Create solver with DG operators
p = 3  # polynomial degree (N needs to be divisible by p + 1)
D_legendre = legendre_derivative_operator(-1.0, 1.0, p + 1)
uniform_mesh = UniformPeriodicMesh1D(coordinates_min, coordinates_max, div(N, p + 1))

D1 = couple_discontinuously(D_legendre, uniform_mesh)
D_pl = couple_discontinuously(D_legendre, uniform_mesh, Val(:plus))
D_min = couple_discontinuously(D_legendre, uniform_mesh, Val(:minus))
D2 = sparse(D_pl) * sparse(D_min)

solver = Solver(D1, D2)
```

This approach creates a DG discretization with polynomial degree `p` within each element. Similarly, you can create a CG discretization using `couple_continuously`. Note that in this case, `N` needs to be divisible by `p` instead of `p + 1`.

## Fourier Spectral Methods

[Fourier collocation methods](@ref fourier_sbp) can be interpreted as periodic SBP operators, which can be constructed via `fourier_derivative_operator`:

```julia
using SummationByPartsOperators: fourier_derivative_operator

# Create solver with Fourier pseudospectral collocation method
D1 = fourier_derivative_operator(xmin(mesh), xmax(mesh), nnodes(mesh))
solver = Solver(D1)
```

## Variable Coefficient Operators

Variable coefficient operators can be useful when solving PDEs that involve elliptic systems with variable coefficients. Using them can lead to fewer allocations.

```julia
using SummationByPartsOperators: Mattsson2012, derivative_operator,
                                 var_coef_derivative_operator

accuracy_order = 2
D1 = derivative_operator(Mattsson2012();
                         derivative_order = 1, accuracy_order,
                         xmin = xmin(mesh), xmax = xmax(mesh), N = N)
# Create a variable-coefficient second-derivative operator
# Initialize with coefficient function `one` - coefficients will be set during simulation
D2 = var_coef_derivative_operator(Mattsson2012(),
                                  2, accuracy_order,
                                  xmin(mesh), xmax(mesh), N, one)
solver = Solver(D1, D2)
```

!!! warning "Limited compatibility"
    Variable coefficient operators (`VarCoefDerivativeOperator`) are currently only supported for the Serre-Green-Naghdi equations with reflecting boundary conditions and flat bathymetry. This is a specialized feature for specific use cases.

## Using Your Custom Solver

Once you have created your solver object using any of the approaches above, you can use it in a `Semidiscretization`:

```julia
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions)
```

The solver you choose should be compatible with your boundary conditions:

- Periodic operators (created with `periodic_derivative_operator` or `fourier_derivative_operator`) require `boundary_condition_periodic`
- Non-periodic operators (created with `MattssonNordström2004`, etc.) are needed for `boundary_condition_reflecting`

## Common Pitfalls

### Mismatched Operators and Boundary Conditions

**Problem**: Using periodic operators with reflecting boundary conditions or vice versa.

```julia
# ❌ This will fail
D1 = periodic_derivative_operator(1, 4, mesh.xmin, mesh.xmax, mesh.N)
solver = Solver(D1)
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_condition_reflecting)  # Error!
```

**Solution**: Match operator type to boundary conditions:

```julia
# ✅ Correct approach
D1 = derivative_operator(MattssonNordström2004(), derivative_order = 1,
                         accuracy_order = 4, xmin = mesh.xmin, xmax = mesh.xmax, N = mesh.N)
solver = Solver(D1)
semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_condition_reflecting)
```

### Incorrect Grid Size for DG Methods

**Problem**: Grid size not divisible by polynomial degree + 1 for DG methods.

```julia
# ❌ N = 101, p = 3, but 101 is not divisible by (3+1) = 4
N = 101
p = 3
mesh = Mesh1D(coordinates_min, coordinates_max, N)  # Error in DG setup!
```

**Solution**: Ensure `N` is divisible by `p + 1`:

```julia
# ✅ Adjust N to be divisible by p + 1
p = 3
N = 100  # 100 ÷ 4 = 25 elements exactly
# Or N = 96, 104, etc. - any multiple of 4
```
