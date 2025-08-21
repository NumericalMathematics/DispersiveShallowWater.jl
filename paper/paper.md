---
title: 'DispersiveShallowWater.jl: A Julia library of structure-preserving numerical methods for dispersive wave equations'
tags:
  - Julia
  - numerical analysis
  - differential equations
authors:
  - name: Joshua Lampert
    orcid: 0009-0007-0971-6709
    affiliation: 1
  - name: Collin Wittenstein
    orcid: 0009-0006-8591-278X
    affiliation: 2
  - name: Hendrik Ranocha
    orcid: 0000-0002-3456-2277
    affiliation: 2
affiliations:
 - name: Department of Mathematics, University of Hamburg, Germany
   index: 1
 - name: Institute of Mathematics, Johannes Gutenberg University Mainz, Germany
   index: 2
date: 21 August 2025
bibliography: paper.bib
---

# Summary

DispersiveShallowWater.jl is a Julia library designed for the numerical simulation of dispersive wave equations, with a focus on structure-preserving methods. The library aims to provide a flexible and efficient
framework for researchers and users working in the field of wave dynamics. The design of DispersiveShallowWater.jl emphasizes broad applicability to a variety of dispersive wave equations, supporting
multiple numerical methods including finite difference (FD), discontinuous Galerkin (DG), continuous Galerkin (CG), and Fourier pseudospectral approaches enabled by the generality of
summation-by-parts (SBP) operators. A central goal is the preservation of key structural properties of the underlying equations, ensuring physically meaningful and robust simulations across diverse scenarios.

# Statement of need

Dispersive wave equations are fundamental in modeling various physical phenomena including shallow water waves, tsunamis, inundations, and other geophysical flows. These phenomena are often modelled as
partial differential equations (PDEs) that exhibit dispersive behavior. The complexity and analytical intractability of most dispersive wave PDEs make the use of numerical methods inevitable for their practical solution.

Accurate simulation of dispersive wave equations requires advanced numerical methods capable of capturing both nonlinear and dispersive effects. It is crucial that these methods preserve fundamental physical
invariants, such as conservation laws and stability, to ensure that the resulting simulations remain physically meaningful, accurate, and reliable. DispersiveShallowWater.jl addresses these needs by implementing
structure-preserving algorithms tailored for a wide range of dispersive wave models.
Summation-by-parts (SBP) operators play a crucial role in the library by enabling the construction of numerical schemes that mimic the integration-by-parts property at the discrete level.
This property is essential for ensuring provable conservation and stability in the numerical solution of dispersive wave equations, making SBP operators particularly well-suited for structure-preserving simulations.
In recent years, SBP operators gained significant attention in the numerical analysis community and have been successfully applied to a variety of problems including dispersive wave equations [@ranocha2021broad; @ranocha2021rate; @linders2023resolving; @lampert2024structure; @ranocha2025structure].

Despite the importance of structure-preserving methods for dispersive wave equations, such approaches are rarely available in existing open-source software packages. DispersiveShallowWater.jl fills this
gap by providing a unified framework that enables users to easily compare different mathematical models, numerical discretizations, and physical setups. This facilitates reproducible research and
accelerates the development and assessment of new methods in the field.

# Features

TODO


# Related research and software

TODO
@lampert2024structure


# Acknowledgements

JL acknowledges the support by the Deutsche Forschungsgemeinschaft (DFG)
within the Research Training Group GRK 2583 "Modeling, Simulation and
Optimization of Fluid Dynamic Applications".
HR TODO


# References
