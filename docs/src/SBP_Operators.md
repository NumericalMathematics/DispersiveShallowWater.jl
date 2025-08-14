# [Summation-by-Parts Operators](@id sbp_operators)

This chapter covers the analytical and mathematical background of summation-by-parts operators in general.

To learn more about different solvers and how to use them in DispersiveShallowWater.jl, go to the chapter about [Solvers](@ref solvers).

## 1. Introduction & Overview

In recent years, summation-by-parts (SBP) operators have gained particular interest in computational mathematics as they allow transferring analytical results from the continuous level to numerical methods in a systematic manner. This is achieved by mimicking integration by parts discretely, which is one of the key ingredients for conservation and stability proofs at the continuous level. In this way, many fundamental analytical properties of hyperbolic-dominated partial differential equations can be obtained in a straightforward manner at the discrete level.[^LampertRanocha2024]

SBP operators were first developed for finite difference methods to mimic stability proofs based on integration by parts as traditionally used in finite element methods. However, exact integration can be impossible or computationally expensive in finite element methods, particularly for complex geometries or nonlinear problems. In this case, SBP formulations can be advantageous since they naturally include a quadrature rule through the mass matrix. In particular, split forms can be used with SBP operators to avoid the need for exact integration while maintaining discrete analogs of important analytical properties such as the chain rule and product rule.

Several classes of numerical methods can be formulated via SBP operators, including finite difference methods, finite volume methods, continuous Galerkin methods, discontinuous Galerkin (DG) methods, and flux reconstruction methods. This unifying framework has made SBP operators a cornerstone of structure-preserving numerical methods across various computational disciplines.

### Why SBP Operators Matter for Dispersive Shallow Water Equations

For dispersive shallow water equations implemented in DispersiveShallowWater.jl, maintaining physical properties like mass and energy conservation is crucial for:

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
    \underbrace{\boldsymbol{u}^T M D_1 \boldsymbol{v}}_{\approx} &+& \underbrace{\boldsymbol{u}^T D_1^T M \boldsymbol{v}}_{\approx} &=& \underbrace{\boldsymbol{u}^T \boldsymbol{e}_R \boldsymbol{e}_R^T \boldsymbol{v}}_{=} &-& \underbrace{\boldsymbol{u}^T \boldsymbol{e}_L \boldsymbol{e}_L^T \boldsymbol{v}}_{=}\\
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

In practice SBP operators come in various *flavours*. 

upwind vs. central vs. FD vs. DG vs. CG vs. Fourier

### Central vs. Upwind SBP Operators

The choice between central and upwind SBP operators is crucial and depends on the nature of your problem.

#### Central SBP Operators


#### [Upwind SBP Operators](@id upwind_sbp)

Upwind operators come in pairs ``D_+`` and ``D_-`` with additional dissipation:
- **Key Property**: ``M(D_+ - D_-)`` is negative semidefinite (provides controlled dissipation)
- **Major Advantage**: **Higher resolution** - better at resolving sharp features
- **Construction**: ``D = (D_+ + D_-)/2`` gives a central operator

#### Why Upwind Operators Have Higher Resolution

The resolution advantage comes from the stencil structure. Consider discretizing a second derivative:

**Central approach**:

```math
(D^2 u)_i \approx \frac{u_{i+2} - 2u_i + u_{i-2}}{4\Delta x^2} \quad \text{[Wide stencil!]}
```

This wide stencil has **grid oscillations** of the form ``(-1,+1,-1,+1,\ldots)`` in its null space.

**Upwind approach**:

```math
(D^2 u)_i \approx \frac{u_{i+1} - 2u_i + u_{i-1}}{\Delta x^2} \quad \text{[Narrow stencil!]}
```

The narrow stencil naturally damps high-frequency oscillations, leading to **better resolution**.

### Different SBP Implementations

#### Finite Difference (FD) SBP

#### [Discontinuous Galerkin (DG) SBP](@id dg_sbp)

#### [Continuous Galerkin (CG) SBP](@id cg_sbp)

#### [Fourier/Spectral SBP](@id fourier_sbp)


