"""
    data_dingemans()

Load the experimental data for the Dingemans experiment. Returns the time values, x-coordinates of the six wave gauges,
and experimental data, which is a matrix of wave heights at each gauge location (columns) over time (rows).

- Dingemans (1994)
  Comparison of computations with Boussinesq-like models and laboratory measurements
  [URL: https://resolver.tudelft.nl/uuid:c2091d53-f455-48af-a84b-ac86680455e9](https://resolver.tudelft.nl/uuid:c2091d53-f455-48af-a84b-ac86680455e9)
- Dingemans (1997):
  Water Wave Propagation Over Uneven Bottoms (In 2 Parts).
  [DOI: 10.1142/1241](https://doi.org/10.1142/1241)

See also: [Dingemans Experiment](https://numericalmathematics.github.io/DispersiveShallowWater.jl/stable/dingemans/) 
in the documentation for more details about this experiment and its use in validation.
"""
function data_dingemans()
    path_dingemans = joinpath(data_dir(), "Dingemans.csv")
    all_data, _ = readdlm(path_dingemans, ','; header = true)
    t_values = all_data[:, 1]
    x_values = (3.04, 9.44, 20.04, 26.04, 30.44, 37.04)
    experimental_data = all_data[:, 2:end]
    return t_values, x_values, experimental_data
end
