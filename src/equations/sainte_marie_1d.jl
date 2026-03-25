@doc raw"""
    SainteMarieEquations1D(; bathymetry_type = bathymetry_mild_slope,
                             gravity,
                             eta0 = 0.0)

Sainte-Marie system in one spatial dimension (with parameter ``\gamma = 2``
compared to the original literature).
The equations are given by
```math
\begin{aligned}
  h_t + (h v)_x &= 0,\\
  h v_t + \frac{1}{2} g (h^2)_x + \frac{1}{2} h (v^2)_x
    + (h p)_x &= -(g h + 2 p) b_x,\\
  h w_t + h v w_x &= 2 p,\\
  v_x + (w - v b_x) / (h / 2) &= 0.
\end{aligned}
```
The unknown quantities of the Sainte-Marie equations are the
total water height ``\eta = h + b`` and the velocity ``v``.
The gravitational acceleration `gravity` is denoted by ``g`` and the bottom topography
(bathymetry) ``b = \eta_0 - D``. The water height above the bathymetry
is therefore given by ``h = \eta - \eta_0 + D``.
The total water height is therefore given by ``\eta = h + b``.
The last two equations of the system determine the auxiliary variable ``w``
and the non-hydrostatic pressure ``p``.

The Sainte-Marie system is similar to the Serre-Green-Naghdi system,
but with a different non-hydrostatic pressure term. In particular, it
also requires solving elliptic equations to determine the time derivative
of the velocity.

Two types of `bathymetry_type` are supported:
- [`bathymetry_flat`](@ref): flat bathymetry (typically ``b = 0`` everywhere)
- [`bathymetry_mild_slope`](@ref): variable bathymetry with mild-slope approximation

!!! note
    The `SainteMarieEquations1D` are only implemented to support their
    hyperbolizaiton [`HyperbolicSainteMarieEquations1D`](@ref) and are not
    intended to be solved directly.

References for the Sainte-Marie system can be found in
- Sainte-Marie (2011)
  Vertically averaged models for the free surface non-hydrostatic
  Euler system: derivation and kinetic interpretation
  [DOI: 10.1142/S0218202511005118](https://doi.org/10.1142/S0218202511005118)
- Bristeau, Mangeney, Sainte-Marie and Seguin (2015)
  An energy-consistent depth-averaged Euler system:
  derivation and properties
  [DOI: 10.3934/dcdsb.2015.20.961](https://doi.org/10.3934/dcdsb.2015.20.961)
- Aïssiouene, Bristeau, Godlewski, Mangeney, Parés Madroñal and Sainte-Marie (2020)
  A two-dimensional method for a family of dispersive shallow water models
  [DOI: 10.5802/smai-jcm.66](https://doi.org/10.5802/smai-jcm.66)
- Escalante, Dumbser and Castro (2019)
  An efficient hyperbolic relaxation system for dispersive non-hydrostatic
  water waves and its solution with high order discontinuous Galerkin schemes
  [DOI: 10.1016/j.jcp.2019.05.035](https://doi.org/10.1016/j.jcp.2019.05.035)
"""
struct SainteMarieEquations1D{Bathymetry <: Union{BathymetryFlat, BathymetryMildSlope},
                              RealT <: Real} <: AbstractSainteMarieEquations{1, 3}
    bathymetry_type::Bathymetry # type of bathymetry
    gravity::RealT # gravitational acceleration
    eta0::RealT # constant still-water surface
end

function SainteMarieEquations1D(; bathymetry_type = bathymetry_mild_slope,
                                gravity,
                                eta0 = 0.0)
    return SainteMarieEquations1D(bathymetry_type, gravity, eta0)
end

# TODO: Consider implementing more functionality to allow comparison
#       with their hyperbolization. For example, a soliton solution
#       is given in https://doi.org/10.3934/dcdsb.2015.20.961.
