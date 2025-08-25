# [Equations](@id equations)

[DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl) provides six different dispersive shallow water equation systems for modeling water waves. Each equation system offers different levels of physical accuracy, computational complexity, and supports various boundary conditions and bathymetry types.

## [Supported Models and Features](@id eq_overview)

The following table provides an overview of all available equation systems and their supported features:

| Equation | Variables | [Periodic boundary conditions](@ref boundary_condition_periodic) | [Reflecting boundary conditions](@ref boundary_condition_reflecting) | [Flat Bathymetry](@ref bathymetry_flat) | [Mild-slope Bathymetry](@ref bathymetry_mild_slope) | [Variable Bathymetry](@ref bathymetry_variable) | Relaxation | Source Terms |
|----------|:---------:|:-----------:|:-------------:|:----:|:-----------:|:--------:|:----------:|:-------:|
| [`KdV`](@ref KdVEquation1D) | ``(\eta)`` | ✅ | ❌ | ✅ | ❌ | ❌ | ✅ | ✅ |
| [`BBM`](@ref BBMEquation1D) | ``(\eta)`` | ✅ | ❌ | ✅ | ❌ | ❌ | ✅ | ✅ |
| [`BBM-BBM`](@ref BBMBBMEquations1D) | ``(\eta, v, D)`` | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ |
| [`Svärd-Kalisch`](@ref SvaerdKalischEquations1D) | ``(\eta, v, D)`` | ✅ | ✅ᵃ | ❌ | ❌ | ✅ | ✅ | ✅ |
| [`Serre-Green-Naghdi`](@ref SerreGreenNaghdiEquations1D) | ``(\eta, v, D)`` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| [`Hyperbolic SGN`](@ref HyperbolicSerreGreenNaghdiEquations1D) |``(\eta, v, D, w, H)`` | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ |

*ᵃReflecting boundary conditions for Svärd-Kalisch equations require `alpha = gamma = 0`*

### Variable Descriptions

![water height and bathymetry](bathymetry.png)

- ``\eta``: Total water height
- ``v``: Velocity in horizontal direction
- ``D``: Still-water depth
- ``w``: Auxiliary variable in hyperbolic approximation (``\approx -h v_x``)
- ``H``: Auxiliary variable in hyperbolic approximation (``\approx h``)

## Abstract Shallow Water Equations Interface

Several equation systems in [DispersiveShallowWater.jl](https://github.com/NumericalMathematics/DispersiveShallowWater.jl) (`BBMBBMEquations1D`, `SvaerdKalischEquations1D`, `SerreGreenNaghdiEquations1D`, and `HyperbolicSerreGreenNaghdiEquations1D`) are subtypes of [`AbstractShallowWaterEquations`](@ref). This design reflects that these systems all contain the classical shallow water equations as a subsystem, extended with additional dispersive terms.

The common interface provides shared functionality like [`waterheight`](@ref), [`velocity`](@ref), [`energy_total`](@ref), and [`entropy`](@ref). This enables consistent analysis and visualization across different dispersive models while maintaining the underlying shallow water physics.

## Hyperbolic Approximations

Some equations are hyperbolic approximations of other systems (e.g., [`HyperbolicSerreGreenNaghdiEquations1D`](@ref) approximates [`SerreGreenNaghdiEquations1D`](@ref)). These systems support two approaches for initial conditions:

1. **Full variables**: Specify all primitive variables including auxiliary variables ``(\eta, v, D, w, H)``
2. **Physical variables**: Specify only the physical variables ``(\eta, v, D)`` from the limit system. The auxiliary variables are automatically initialized using appropriate approximations.

This flexibility allows using the same initial conditions for both the original and hyperbolic systems, simplifying comparisons and testing.

### Detailed Documentation

Each equation system below includes its complete mathematical formulation, physical background, implementation details, conservation properties, and relevant literature references.

## Korteweg–De Vries (KdV)

```@docs
DispersiveShallowWater.KdVEquation1D
```

## Benjamin-Bona-Mahony (BBM)

```@docs
DispersiveShallowWater.BBMEquation1D
```

## BBM-BBM

```@docs
DispersiveShallowWater.BBMBBMEquations1D
```

## Svärd-Kalisch

```@docs
DispersiveShallowWater.SvaerdKalischEquations1D
DispersiveShallowWater.SvärdKalischEquations1D
```

## Serre-Green-Naghdi

```@docs
DispersiveShallowWater.SerreGreenNaghdiEquations1D
```

## Hyperbolic Serre-Green-Naghdi

```@docs
DispersiveShallowWater.HyperbolicSerreGreenNaghdiEquations1D
```
