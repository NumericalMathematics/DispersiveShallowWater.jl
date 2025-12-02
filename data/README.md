# Experimental Data

This directory contains experimental datasets used in common benchmark problems in dispersive wave modeling. Each dataset has a corresponding function that returns the data as Julia objects for easy integration with DispersiveShallowWater.jl simulations.

## Available Datasets

### Dingemans Experiment Data (`Dingemans.csv`)

The Dingemans experiment is a classic benchmark for validating dispersive shallow water models. It features waves propagating over a trapezoidal bathymetry profile, providing experimental data for wave transformation over varying bottom topography.

The data are taken from the experiments of Maarten Dingemans:

```bibtex
@techreport{dingemans1994comparison,
  title={Comparison of computations with {B}oussinesq-like models and laboratory measurements},
  author={Maarten W. Dingemans},
  institution={Delft Hydraulics},
  year={1994},
  number={H1684.12},
  url={http://resolver.tudelft.nl/uuid:c2091d53-f455-48af-a84b-ac86680455e9}
}

@book{dingemans1997water,
  title={Water Wave Propagation Over Uneven Bottoms},
  author={Maarten W. Dingemans},
  year={1997},
  volume={13},
  doi={10.1142/1241},
  publisher={World Scientific}
}
```


#### Data Format

The CSV file contains time series data with the first column being time values and subsequent columns representing wave height measurements at different spatial locations (six wave gauges).

#### Usage

In order to access the data, execute:

```julia
t_values, x_values, experimental_data = data_dingemans()
```

The function returns:
- `t_values`: Time values at which measurements were recorded
- `x_values`: Spatial coordinates of the six wave gauges (3.04, 9.44, 20.04, 26.04, 30.44, 37.04)
- `experimental_data`: Matrix of wave heights at each gauge location (columns) over time (rows)

For more details about this experiment and its use in validation, see the [Dingemans Experiment](https://numericalmathematics.github.io/DispersiveShallowWater.jl/stable/dingemans/) section in the documentation.