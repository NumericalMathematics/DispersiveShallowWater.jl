# [Summation-by-Parts Operators](@id sbp_operators)

This chapter covers the analytical and mathematical background of summation-by-parts operators in general.

To learn more about different solvers and how to use them in [DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl), go to the chapter about [Solvers](@ref solvers).

## 1. Introduction & Overview

In recent years, summation-by-parts (SBP) operators have gained particular interest in computational mathematics as they allow transferring analytical results from the continuous level to numerical methods in a systematic manner. This is achieved by mimicking integration by parts discretely, which is one of the key ingredients for conservation and stability proofs at the continuous level. In this way, many fundamental analytical properties of hyperbolic-dominated partial differential equations can be obtained in a straightforward manner at the discrete level.[^LampertRanocha2024]

SBP operators were first developed for finite difference methods to mimic stability proofs based on integration by parts as traditionally used in finite element methods. However, exact integration can be impossible or computationally expensive in finite element methods, particularly for complex geometries or nonlinear problems. In this case, SBP formulations can be advantageous since they naturally include a quadrature rule through the mass matrix. In particular, split forms can be used with SBP operators to avoid the need for exact integration while maintaining discrete analogs of important analytical properties such as the chain rule and product rule.

Several classes of numerical methods can be formulated via SBP operators, including finite difference methods, finite volume methods, continuous Galerkin methods, discontinuous Galerkin (DG) methods, and flux reconstruction methods. This unifying framework has made SBP operators a cornerstone of structure-preserving numerical methods across various computational disciplines.

### Why SBP Operators Matter for Dispersive Shallow Water Equations

For dispersive shallow water equations implemented in [DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl), maintaining physical properties like mass and energy conservation is crucial for:

- **Long-time stability**: Simulations remain stable over extended time periods without spurious growth of numerical errors
- **Physical accuracy**: The discrete solution respects fundamental physical laws such as conservation of mass, momentum, and energy
- **Robustness**: Methods are less sensitive to parameter choices and grid resolution, leading to more reliable simulations

SBP operators achieve this by providing:

- **Exact conservation**: Discrete conservation laws that hold to machine precision
- **Provable stability**: Mathematical guarantees about the behavior of the numerical method through discrete energy estimates
- **Flexibility**: A unified framework that encompasses finite differences, finite elements, and spectral methods

[^LampertRanocha2024]:
    Lampert, Ranocha (2024):
    Structure-Preserving Numerical Methods for Two Nonlinear Systems of Dispersive Wave Equations
    [arXiv: 2402.16669](https://arxiv.org/abs/2402.16669)

## 2. Mathematical Foundation

### The Core SBP Property

The fundamental property that defines a first-derivative SBP operator is:

```math
MD + D^T M = \boldsymbol{e}_R \boldsymbol{e}_R^T - \boldsymbol{e}_L \boldsymbol{e}_L^T
```

where:

- ``D`` is the derivative operator (matrix)
- ``M`` is the symmetric, positive definite mass matrix
- ``\boldsymbol{e}_L = (1,0,...,0)^T`` and ``\boldsymbol{e}_R = (0,...,0,1)^T`` extract boundary values

This property is the discrete analog of integration by parts:

```math
\begin{array}{ccccccc}
    \underbrace{\boldsymbol{u}^T M D \boldsymbol{v}}_{\approx} &+& \underbrace{\boldsymbol{u}^T D^T M \boldsymbol{v}}_{\approx} &=& \underbrace{\boldsymbol{u}^T \boldsymbol{e}_R \boldsymbol{e}_R^T \boldsymbol{v}}_{=} &-& \underbrace{\boldsymbol{u}^T \boldsymbol{e}_L \boldsymbol{e}_L^T \boldsymbol{v}}_{=}\\
    \overbrace{\displaystyle\int_\Omega u v_x \, dx} &+& \overbrace{\displaystyle\int_\Omega u_x v \, dx} &=& \overbrace{u(x_{\text{max}}) v(x_{\text{max}})} &-& \overbrace{u(x_{\text{min}}) v(x_{\text{min}})}.
\end{array}
```

where for periodic SBP operators their property naturally simplifies to

```math
MD + D^T M = 0.
```

### Understanding the Mass Matrix

The mass matrix ``M`` approximates the continuous inner product:

```math
\langle \boldsymbol{u}, \boldsymbol{v} \rangle_M = \boldsymbol{u}^T M \boldsymbol{v} \approx \int_{x_{min}}^{x_{max}} u(x) v(x) dx
```

For the approximation to be meaningful, we require:

```math
\boldsymbol{1}^T M \boldsymbol{1} = x_{\max} - x_{\min}
```

This ensures that the discrete inner product correctly integrates constants. This property is always fulfilled for at least first-order SBP operators.

### Understanding the Derivative Operator

The derivative operator ``D`` approximates the spatial derivative of a function on a discrete grid. That is, for a smooth function ``u(x)``, we want the discrete expression ``D \boldsymbol{u}`` to approximate ``\partial_x u(x)`` at the grid points.

A first-derivative SBP operator is said to be **$p$-th order accurate** if it satisfies:

```math
D\, \boldsymbol{x}^k = k\, \boldsymbol{x}^{k-1}, \quad \text{for all } k = 0, 1, ..., p,
```

where ``\boldsymbol{x}^k = (x_1^k, \dots, x_N^k)^T``. This means the discrete operator correctly differentiates monomials up to degree $p$. In particular, **consistency** requires at least zeroth-order accuracy, i.e.,

```math
D\, \boldsymbol{1} = \boldsymbol{0}.
```

However, approximating derivatives near boundaries poses a challenge: interior points can use standard finite difference stencils, but near the edges (e.g., at `x_1` or `x_N`), one must use specially designed one-sided approximations that still preserve stability and accuracy. Periodic SBP Operators do not have this problem.

### Example: Central Difference Operator

A simple example is the classical second-order central difference operator on a uniform grid:

```math
D = \frac{1}{\Delta x} \begin{pmatrix}
-1 & 1 & & & \\
-1/2 & 0 & 1/2 & & \\
& -1/2 & 0 & 1/2 & \\
& & \ddots & \ddots & \ddots \\
& & & -1 & 1
\end{pmatrix}
```

with mass matrix:

```math
M = \Delta x \cdot \text{diag}(1/2, 1, 1, \ldots, 1, 1/2)
```

You can verify that this satisfies the SBP property and provides second-order accuracy in the interior with first-order accuracy at the boundaries.

## 3. Types of SBP Operators

In practice SBP operators come in various *flavors*. Consult the [documentation of SummationByPartsOperators.jl](https://ranocha.de/SummationByPartsOperators.jl/stable/)
for more details on how to construct the different types.

### Central Finite Difference (FD) SBP Operators

The simplest form of SBP operators are central finite difference operators and these are the operators that are mostly used in [DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl)
For periodic boundary conditions, they can be created with `periodic_derivative_operator` and for non-periodic boundary conditions with `derivative_operator`
from [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl/).

### [Upwind SBP Operators](@id upwind_sbp)

Upwind operators come in pairs ``D_+`` and ``D_-``, which satisfy the relation

```math
MD_+ + D_-^T M = \boldsymbol{e}_R \boldsymbol{e}_R^T - \boldsymbol{e}_L \boldsymbol{e}_L^T,
```

where additionally ``M(D_+ - D_-)`` is negative semidefinite. This property is often useful to construct dissipative numerical schemes, which have an improved stability.
The operators ``D_+`` and ``D_-`` are biased in one direction and are therefore useful in simulations with unidirectional flow, which favor a specific flow direction.
They can also be helpful to construct (central) second-derivative operators by ``D_2 = D_+D_-`` (or ``D_2 = D_-D_+``), which is often advantageous compared to wide-stencil
central FD operators like ``D_2 = D_1^2`` because second-derivative operators based on upwind operators have a narrow stencil leading to a better resolution.
With [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl/), upwind operators can be constructed using `upwind_operators`. Each part of the operator can be accessed by the `minus` and `plus` fields.
Additionally, upwind operators induce a central first-derivative operator by ``D = (D_+ + D_-)/2``, which can be accessed by `central`.

### [Discontinuous Galerkin (DG) SBP](@id dg_sbp)

Classically, SBP operators were developed from the perspective of finite difference methods. However, more recently especially in [^Gassner2013] and subsequent papers,
the connection of finite element method to SBP operators was developed. It turns out that many finite element schemes, like the discontinuous Galerkin spectral element
method (DGSEM), can be interpreted as schemes based on certain SBP operators. In this context, the SBP operators are defined locally on a reference element. The
DGSEM uses Gauss-Lobatto-Legendre nodes and weights to form a quadrature rule, which also leads to a differentiation matrix. Together, these satisfy the SBP property.
To obtain global mass and derivative matrices, the local operators can be coupled across the elements as presented in [^RanochaMitsotakisKetcheson2021]. You can create periodic
DG operators with

```julia
coordinates_min, coordinates_max = -1.0, 1.0
N = 16 # N needs to be divisible by p + 1, i.e. this corresponds to N/(p + 1) = 4 elements
p = 3 # polynomial degree
D_legendre = legendre_derivative_operator(-1.0, 1.0, p + 1)
uniform_mesh = UniformPeriodicMesh1D(coordinates_min, coordinates_max, div(N, p + 1))
D = couple_discontinuously(D_legendre, uniform_mesh)
```

Using `couple_discontinuously`, you can also construct upwind SBP operators by additionally passing `Val(:plus)` or `Val(:minus)`. Note that this construction results into non-uniformly
distributed nodes and due to the discontinuous nature to repeated nodes at interfaces between elements.

### [Continuous Galerkin (CG) SBP](@id cg_sbp)

Similarly to DG operators, CG operators can be constructed using `couple_continuously`. In contrast to the DG operators, `N` needs to be divisible by `p` and there are no repeated nodes at the interfaces.

### [Fourier/Spectral SBP](@id fourier_sbp)

Fourier or spectral SBP operators are constructed using Fourier basis functions. These operators can be used for problems with periodic boundary conditions. The key idea is to represent the solution in
terms of its Fourier coefficients and to apply differentiation in the Fourier space. Note that these operators have dense derivative matrices and are therefore often more computationally expensive. In
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl/), they can be constructed with `fourier_derivative_matrix`.

### Variable Coefficient Operators

A special class of SBP operators is given by variable coefficient operators, which are discrete operators for the second derivative approximating terms of the form ``\partial_x(b \partial_x u)``. Directly
incorporating the variable coefficient `b` into the SBP operator is desirable compared to subsequent application of first-derivative operators ``D \textrm{diag}(\boldsymbol{b}) D`` because it leads to a
more compact stencil and therefore improved numerical properties. You can use `var_coef_derivative_operator` with source `Mattsson2012` to construct such operators.

[^Gassner2013]:
    Gassner (2013):
    A Skew-Symmetric Discontinuous Galerkin Spectral Element Discretization
    and Its Relation to SBP-SAT Finite Difference Methods
    [DOI: 10.1137/120890144](https://epubs.siam.org/doi/10.1137/120890144)

[^RanochaMitsotakisKetcheson2021]:
    Ranocha, Mitsokatis, Ketcheson (2021):
    A broad class of conservative numerical methods for dispersive wave equations
    [DOI: 10.4208/cicp.oa-2020-0119](https://doi.org/10.4208/cicp.oa-2020-0119)
