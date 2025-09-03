This directory collects some experimental data used in common benchmark problems in dispersive wave modeling. Each set of data has a corresponding function, which returns the data as Julia objects.

# Dingemans.csv


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

The positions x1 to x6 are (3.04, 9.44, 20.04, 26.04, 30.44, 37.04) respectively. See `data_dingemans`.