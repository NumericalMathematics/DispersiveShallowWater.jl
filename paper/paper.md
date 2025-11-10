---
title: 'DispersiveShallowWater.jl: A Julia library of structure-preserving numerical methods for dispersive wave equations'
tags:
  - Julia
  - numerical analysis
  - differential equations
  - dispersive wave equations
  - summation-by-parts
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

Dispersive wave equations are fundamental in modeling various physical phenomena including shallow water waves, tsunamis, inundations, and other geophysical flows. These phenomena are often modeled as
partial differential equations (PDEs) that exhibit dispersive behavior. The complexity and analytical intractability of most dispersive wave PDEs make the use of numerical methods inevitable for their practical solution.

Accurate simulation of dispersive wave equations requires advanced numerical methods capable of capturing both nonlinear and dispersive effects. It is crucial that these methods preserve fundamental physical
invariants, such as conservation laws and stability, to ensure that the resulting simulations remain physically meaningful, accurate, and reliable. DispersiveShallowWater.jl addresses these needs by implementing
structure-preserving algorithms tailored for a wide range of dispersive wave models.
Summation-by-parts (SBP) operators play a crucial role in the library by enabling the construction of numerical schemes that mimic the integration-by-parts property at the discrete level.
This property is essential for ensuring provable conservation and stability in the numerical solution of dispersive wave equations, making SBP operators particularly well-suited for structure-preserving simulations.
In recent years, SBP operators have gained significant attention in the numerical analysis community and have been successfully applied to a variety of problems including dispersive wave equations [@mattsson2014diagonal; @almquist2014high; @mattsson2016high; @lindeberg2021high; @ranocha2021broad; @ranocha2021rate; @rydin2021high; @linders2023resolving; @lampert2024structure; @biswas2025traveling; @kjelldahl2025numerical; @ranocha2025structure; @giesselmann2025convergence].

Despite the importance of structure-preserving methods for dispersive wave equations, such approaches are rarely available in existing open-source software packages. DispersiveShallowWater.jl is
specifically designed to serve researchers who develop and compare numerical algorithms and mathematical models for dispersive wave phenomena. By providing a unified framework, the library enables
users to systematically evaluate different models, numerical discretizations, and physical setups. This focus facilitates reproducible research and
accelerates the development and assessment of new methods in the field.

# Features

DispersiveShallowWater.jl is written in the Julia programming language [@bezanson2017julia] and leverages Julia's strengths in scientific computing, such as high performance, ease of use,
and rich ecosystem of libraries and tools.

To date, DispersiveShallowWater.jl supports classical one-dimensional scalar dispersive wave equations like the Korteweg-de Vries (KdV) equation [@korteweg1895change] and the
Benjamin-Bona-Mahony (BBM) equation [@benjamin1972model] as well as more sophisticated one-dimensional systems of equations like the BBM-BBM system [@bona1998boussinesq], the Serre-Green-Naghdi equations
[@serre1953contribution; @green1976derivation], a hyperbolic approximation thereof [@favrie2017rapid], and the Svärd-Kalisch equations [@svard2025novel].

The package integrates well into the existing ecosystem of Julia using SummationByPartsOperators.jl [@ranocha2021sbp] for the construction of SBP operators,
OrdinaryDiffEq.jl [@rackauckas2017differentialequations] for solving the ordinary differential equations resulting from spatial discretization, and Plots.jl [@christ2023plots] for visualization. This allows the library to use advanced techniques
implemented in these packages and benefit from their extensive functionality. Moreover, the design concept behind DispersiveShallowWater.jl is largely inspired by the well-established numerical
solution framework Trixi.jl [@ranocha2022adaptive; @schlottkelakemper2021purely; @schlottkelakemper2025trixi] making the interface familiar to users of Trixi.jl and easy to extend.

In addition, users benefit from a suite of built-in analysis and postprocessing tools for investigation of numerical and physical properties, performance evaluation, and visualization. Entropy-conserving
time integration schemes based on relaxation approaches are implemented, enabling stability also on the fully-discrete level, see @ketcheson2019relaxation, @ranocha2020relaxation.
Furthermore, DispersiveShallowWater.jl includes routines for computing and analyzing linear dispersion relations, enabling theoretical investigation and comparison of different models.

# Related research and software

Over the last century, several mathematical models describing the behavior of water waves have been proposed. As, e.g., outlined by @glimsdal2013dispersion, the ability to model dispersion effects is essential
for many applications in fluid dynamics, coastal engineering, and environmental science. Therefore, many equations have been developed to capture these effects, which makes them physically more accurate
compared to, e.g., the well-known shallow water equations, but also numerically more challenging to solve. Hence, researchers have developed a wide range of numerical methods to solve these equations,
including finite difference, finite volume, discontinuous Galerkin, and spectral methods. However, many of these methods do not preserve the underlying structure of the equations, which can lead to
numerical artifacts and inaccuracies in the simulations.

This leads to the recent trend in numerical analysis to develop structure-preserving discretization methods that maintain the physical properties of the equations.
For the dispersive wave equations mentioned above, several structure-preserving methods have been proposed and analyzed in the literature, including the use of summation-by-parts (SBP) operators,
cf. @biswas2025traveling, @ranocha2021broad, @linders2023resolving, @lampert2024structure, and @ranocha2025structure. DispersiveShallowWater.jl provides a unified framework, which offers access
to the numerical discretizations developed in these works. In @lampert2024structure, DispersiveShallowWater.jl is used for the implementation of the presented methods. The work also compares the numerical
solutions to data obtained from experiments showing good agreement.

To the authors' knowledge, no other software package provides the same level of functionality for simulating dispersive shallow water waves as DispersiveShallowWater.jl.
Other open-source software packages, such as
Basilisk ([http://basilisk.fr/](http://basilisk.fr/), accessed 2025-08-22),
FUNWAVE-TVD ([https://fengyanshi.github.io](https://fengyanshi.github.io), accessed 2025-08-22),
SWASH ([https://www.tudelft.nl/swash](https://www.tudelft.nl/swash), accessed 2025-08-22),
and the HySEA family of codes ([https://edanya.uma.es/hysea/](https://edanya.uma.es/hysea/), accessed 2025-08-22)
focus on operational, production-grade simulations while DispersiveShallowWater.jl is aimed at developing and comparing models and numerical methods.
The widely used coastal/tsunami models mentioned above offer similar capabilities for some of the relevant equations but may not include all the features and tools available in DispersiveShallowWater.jl and rely on different numerical methods and approaches.
Other Julia packages, such as Oceananigans.jl [@ramadhan2020oceananigans], GeophysicalFlows.jl [@GeophysicalFlowsJOSS],
SpectralWaves.jl [@paprota2025spectralwaves], WaterWaves1D.jl [@duchenewaterwaves1d], ShallowWaters.jl [@klowershallowwaters], and TrixiShallowWater.jl [@winters2025trixi],
also provide tools for simulating shallow water and dispersive wave phenomena, but differ in their focus, supported models, and numerical methods.
While some research papers offer supplementary code, these are typically limited to small scripts intended for reproducing specific results and are not
developed as general-purpose software libraries.

# Acknowledgments

JL acknowledges the support by the Deutsche Forschungsgemeinschaft (DFG)
within the Research Training Group GRK 2583 "Modeling, Simulation and
Optimization of Fluid Dynamic Applications".
HR additionally acknowledges support from the DFG through individual research grants 513301895 and 528753982, as well as within the DFG priority program SPP 2410 with project number 526031774.


# References
