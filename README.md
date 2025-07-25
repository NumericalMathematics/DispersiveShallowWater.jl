# DispersiveShallowWater.jl

[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://NumericalMathematics.github.io/DispersiveShallowWater.jl/stable)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://NumericalMathematics.github.io/DispersiveShallowWater.jl/dev/)
[![Build Status](https://github.com/NumericalMathematics/DispersiveShallowWater.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/NumericalMathematics/DispersiveShallowWater.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/NumericalMathematics/DispersiveShallowWater.jl/graph/badge.svg)](https://codecov.io/gh/NumericalMathematics/DispersiveShallowWater.jl)
[![Coveralls](https://coveralls.io/repos/github/NumericalMathematics/DispersiveShallowWater.jl/badge.svg?branch=main)](https://coveralls.io/github/NumericalMathematics/DispersiveShallowWater.jl?branch=main)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/635090135.svg)](https://zenodo.org/doi/10.5281/zenodo.10034636)

**DispersiveShallowWater.jl** is a [Julia](https://julialang.org/) package that implements structure-preserving numerical methods for dispersive shallow water models.
To date, it provides provably conservative, entropy-conserving and well-balanced numerical schemes for some dispersive shallow water models:

* the [Korteweg–De Vries (KdV) equation as the prototypical example of an integrable PDE](https://doi.org/10.1007/s10915-025-02898-x),
* the [Benjamin-Bona-Mahony (BBM) equation, also known as regularized long-wave equation](https://doi.org/10.4208/cicp.OA-2020-0119),
* the [BBM-BBM equations with varying bottom topography](https://iopscience.iop.org/article/10.1088/1361-6544/ac3c29),
* the [dispersive shallow water model proposed by Magnus Svärd and Henrik Kalisch](https://arxiv.org/abs/2302.09924),
* the [Serre-Green-Naghdi equations in standard and hyperbolic form](https://arxiv.org/abs/2408.02665).

The semidiscretizations are based on summation-by-parts (SBP) operators, which are implemented in [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl/).
To obtain fully discrete schemes, the time integration methods from [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) are used to solve the resulting ordinary differential equations.
Fully discrete entropy-conservative methods can be obtained by using the [relaxation method](https://epubs.siam.org/doi/10.1137/19M1263662) provided by DispersiveShallowWater.jl.
A more detailed documentation can be found [online](https://NumericalMathematics.github.io./DispersiveShallowWater.jl/stable/).

## Installation

If you have not yet installed Julia, then you first need to [download Julia](https://julialang.org/downloads/). Please [follow the instructions for your operating system](https://julialang.org/downloads/platform/).
DispersiveShallowWater.jl works with Julia v1.10 and newer. DispersiveShallowWater.jl is a registered Julia package. Therefore, you can install it by executing the following commands from the Julia REPL

```julia
julia> using Pkg

julia> Pkg.add(["DispersiveShallowWater", "OrdinaryDiffEqTsit5", "Plots"])
```

In addition, this installs the packages OrdinaryDiffEqTsit5.jl from [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl)
used for time-integration and [Plots.jl](https://github.com/JuliaPlots/Plots.jl) to visualize the results. If you want to use
other time integration methods than `Tsit5`, you can install the respective subpackage or OrdinaryDiffEq.jl, which will install
every available solver.
If you want to use other SBP operators than the default operators that DispersiveShallowWater.jl uses, then you also need [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl),
which can be installed running

```julia
julia> Pkg.add("SummationByPartsOperators")
```

## Usage

In the Julia REPL, first load the package DispersiveShallowWater.jl

```julia
julia> using DispersiveShallowWater
```

You can run a basic simulation that solves the BBM-BBM equations by executing

```julia
julia> include(default_example());
```

The result can be visualized by using the package Plots.jl

```julia
julia> using Plots
julia> plot(semi => sol)
```

The command `plot` expects a `Pair` consisting of a `Semidiscretization` and an `ODESolution`. The visualization can also be customized, see the [documentation](https://NumericalMathematics.github.io/DispersiveShallowWater.jl/stable/overview#visualize_results)
for more details. Other examples can be found in the subdirectory [examples/](https://github.com/NumericalMathematics/DispersiveShallowWater.jl/tree/main/examples).
A list of all examples is returned by running `get_examples()`. You can pass the filename of one of the examples or your own simulation file to `include` in order to run it,
e.g., `include(joinpath(examples_dir(), "svaerd_kalisch_1d", "svaerd_kalisch_1d_dingemans_relaxation.jl"))`.

## Referencing

You can directly refer to DispersiveShallowWater.jl as

```bibtex
@misc{lampert2025dispersive,
  title={{D}ispersive{S}hallow{W}ater.jl: {S}tructure-preserving numerical
         methods for dispersive shallow water models},
  author={Lampert, Joshua and Wittenstein, Collin and Ranocha, Hendrik},
  year={2025},
  howpublished={\url{https://github.com/NumericalMathematics/DispersiveShallowWater.jl}},
  doi={10.5281/zenodo.10034636}
}
```

## Authors

The package is mainly developed and maintained by Joshua Lampert (University of Hamburg)
with contributions from Hendrik Ranocha (Johannes Gutenberg University Mainz)
and Collin Wittenstein (Johannes Gutenberg University Mainz).
Some parts of this repository are based on parts of
[Dispersive-wave-schemes-notebooks. A Broad Class of Conservative Numerical Methods for Dispersive Wave Equations](https://github.com/ranocha/Dispersive-wave-schemes-notebooks)
by Hendrik Ranocha, Dimitrios Mitsotakis and David Ketcheson.
The code structure is inspired by [Trixi.jl](https://github.com/trixi-framework/Trixi.jl/).

## License and contributing

DispersiveShallowWater.jl is published under the MIT license (see [License](https://github.com/NumericalMathematics/DispersiveShallowWater.jl/blob/main/LICENSE)). We are pleased to accept contributions from everyone, preferably in the form of a PR.
