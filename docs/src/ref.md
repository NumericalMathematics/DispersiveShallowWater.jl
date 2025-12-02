# DispersiveShallowWater.jl API

```@meta
CurrentModule = DispersiveShallowWater
```

```@autodocs
Modules = [DispersiveShallowWater]
Pages = ["DispersiveShallowWater.jl"]
```

## Model specific equations

```@autodocs
Modules = [DispersiveShallowWater]
Pages = Main.EQUATIONS_FILES
Filter = t -> !(t in [DispersiveShallowWater.KdVEquation1D,
                     DispersiveShallowWater.BBMEquation1D,
                     DispersiveShallowWater.BBMBBMEquations1D,
                     DispersiveShallowWater.SvaerdKalischEquations1D,
                     DispersiveShallowWater.SerreGreenNaghdiEquations1D,
                     DispersiveShallowWater.HyperbolicSerreGreenNaghdiEquations1D])
```

## Linear dispersion relations

```@autodocs
Modules = [DispersiveShallowWater]
Pages = ["dispersion_relation.jl"]
```

## Mesh

```@autodocs
Modules = [DispersiveShallowWater]
Pages = ["mesh.jl"]
```

## Boundary conditions

```@autodocs
Modules = [DispersiveShallowWater]
Pages = ["boundary_conditions.jl"]
```

## Solver

```@autodocs
Modules = [DispersiveShallowWater]
Pages = ["solver.jl"]
```

## Semidiscretization

```@autodocs
Modules = [DispersiveShallowWater]
Pages = ["semidiscretization.jl"]
```

## Callbacks

```@autodocs
Modules = [DispersiveShallowWater]
Pages = Main.CALLBACKS_STEP_FILES
```

## Utilities

### Experimental data

```@autodocs
Modules = [DispersiveShallowWater]
Pages = ["experimental_data.jl"]
```

### Utilities from DispersiveShallowWater.jl

```@autodocs
Modules = [DispersiveShallowWater]
Pages = ["util.jl"]
```

### Utilities from TrixiBase.jl

Be aware that only `trixi_include` is being exported from DispersiveShallowWater.jl. To access the other [TrixiBase.jl](https://github.com/trixi-framework/TrixiBase.jl) functions, you need to either:

- Use the fully qualified name: `DispersiveShallowWater.timer()`, `DispersiveShallowWater.@trixi_timeit`, etc.
- Import [TrixiBase.jl](https://github.com/trixi-framework/TrixiBase.jl) explicitly: `using TrixiBase` or `import TrixiBase`

```@meta
CurrentModule = TrixiBase
```

```@autodocs
Modules = [TrixiBase]
```
