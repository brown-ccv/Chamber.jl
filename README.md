[![Build Status](https://github.com/brown-ccv/Chamber.jl/actions/workflows/test.yml/badge.svg)](https://github.com/brown-ccv/Chamber.jl/actions?query=workflows/test)
[![][docs-stable-img]][docs-stable-url]
[![][docs-dev-img]][docs-dev-url]

# Chamber.jl
`Chamber.jl` is a Julia package for simulating the eruption of a volcano using a model for the frequency of eruptions of upper crustal magma chambers based on Degruyter and Huber (2014). The package generates a CSV file and figures with the solution data for various variables over time.

## Installation
To install the `Chamber` package, open Julia and use the package manager:

```julia
using Pkg
Pkg.add("Chamber")
```

## Usage
The main function in the `Chamber` package is `chamber`, which simulates the evolution of a magma chamber over time and returns a `DataFrame` with the solution data. The function takes the following arguments:
```julia
chamber(composition, end_time, log_volume_km3, InitialConc_H2O, InitialConc_CO2, log_vfr, depth, output_dirname; kwargs...)
```
### Arguments
- `composition`: The magma composition. Use `Silicic()` for rhyolite composition (arc magma) or `Mafic()` for basalt composition (rift).
- `end_time`: Maximum magma chamber evolution duration in seconds.
- `log_volume_km3`: The initial volume of the chamber in logarithmic scale. The actual initial chamber volume is calculated as 10^(log_volume_km3) in km³.
- `InitialConc_H2O`: The initial weight fraction of water in the magma (exsolved + dissolved).
- `InitialConc_CO2`: The initial weight fraction of CO₂ in the magma (exsolved + dissolved).
- `log_vfr`: Magma recharge rate in km³/yr calculated as 10^(`log_vfr`).
- `depth`: Depth of the magma chamber in meters.
- `output_dirname`: (Optional) Name of the output directory. Defaults to current timestamp.

### Keyword Arguments
- `plotfig`: (Optional, default: `true`). Generate and plot figures for each result if true.

### Returns
A `DataFrame` containing the solution with columns:
- `time`: Simulation timestamps (sec).
- `P+dP`: Pressure (Pa).
- `T`: Temperature (K).
- `eps_g`: Gas volume fraction.
- `V`: Volume of the magma chamber (m³).
- `rho_m`: Density of the melt (kg/m³).
- `rho_x`: Density of magma crystal (kg/m³).
- `X_CO2`: Mole fraction of CO2 in the gas.
- `total_mass`: Total mass of magma chamber (kg).
- `total_mass_H2O`: Total mass of water in the magma (kg).
- `total_mass_CO2`: Total mass of CO₂ in the magma (kg).
- `eps_x`: Crystal volume fraction.

### Outputs
A directory named after `output_dirname` or the default value, containing the following files:
- `out.csv`: A comma-delimited ascii file containing the solution columns listed above.
- `eruptions.csv`: A comma-delimited ascii file. It is organized as follows:
  - `time_of_eruption` (sec)
  - `duration_of_eruption` (sec)
  - `mass_erupted` (kg)
  - `volume_erupted` (km³).
- Figures for P+dP(t), T(t), eps_g(t), V(t), X_CO2(t), total_mass(t).

## Examples
### Example 1: Rhyolite composition (arc magma)
Run a simulation with silicic magma chamber:
```julia
julia> using Chamber

julia> # Define simulation parameters
julia> composition = Silicic()
julia> end_time = 3e9
julia> log_volume_km3 = 0.2
julia> InitialConc_H2O = 0.04
julia> InitialConc_CO2 = 0.001
julia> log_vfr = -3.3
julia> depth = 8e3

julia> # Run simulation
julia> dataframe = chamber(composition, end_time, log_volume_km3, InitialConc_H2O, InitialConc_CO2, log_vfr, depth)
470×12 DataFrame
 Row │ time            P+dP       T        eps_g       V          rho_m    rho_x    X_CO2     total_mass  total_mass_H2O  total_mass_CO2  eps_x    
     │ Float64         Float64    Float64  Float64     Float64    Float64  Float64  Float64   Float64     Float64         Float64         Float64  
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │      0.0        2.1582e8   1046.71  0.00558159  1.58489e9  2400.0   2600.0   0.346043  3.83394e12      1.54039e11       3.85098e9  0.146696
   2 │ 100000.0        2.15824e8  1046.71  0.00558131  1.58489e9  2400.0   2600.0   0.34605   3.83395e12      1.54039e11       3.85099e9  0.146696
   3 │      7.61993e5  2.15848e8  1046.71  0.0055795   1.5849e9   2400.01  2600.01  0.346101  3.83397e12      1.5404e11        3.85101e9  0.146698
  ⋮  │       ⋮             ⋮         ⋮         ⋮           ⋮         ⋮        ⋮        ⋮          ⋮             ⋮               ⋮            ⋮
 468 │      2.98643e9  2.23399e8  1043.44  0.00559005  1.58782e9  2401.9   2602.06  0.33863   3.85078e12      1.54716e11       3.8679e9   0.167628
 469 │      2.99643e9  2.23759e8  1043.43  0.00556433  1.58788e9  2401.99  2602.15  0.339346  3.85116e12      1.54731e11       3.86828e9  0.16766
 470 │      3.0e9      2.23888e8  1043.43  0.00555519  1.5879e9   2402.02  2602.18  0.339601  3.8513e12       1.54737e11       3.86842e9  0.167671
                                                                                                                                   464 rows omitted
```
As noted, the `chamber` function returns a DataFrame containing the solution data with columns described above. Additionally, it automatically creates a directory named with the current timestamp (by default) to store the output files including figures and CSV data files.

### Example 2: Basalt composition (rift magma)
Run a simulation with mafic magma chamber, with custom directory name "MyDirname":
```julia
julia> composition = Mafic()
julia> end_time = 3e9
julia> log_volume_km3 = 0.2
julia> InitialConc_H2O = 0.01
julia> InitialConc_CO2 = 0.001
julia> log_vfr = -3.3
julia> depth = 8e3
julia> output_dirname = "MyDirname"

julia> chamber(composition, end_time, log_volume_km3, InitialConc_H2O, InitialConc_CO2, log_vfr, depth, output_dirname)
458×12 DataFrame
 Row │ time            P+dP       T        eps_g       V          rho_m    rho_x    X_CO2     total_mass  total_mass_H2O  total_mass_CO2  eps_x    
     │ Float64         Float64    Float64  Float64     Float64    Float64  Float64  Float64   Float64     Float64         Float64         Float64  
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │      0.0        2.1582e8   1405.74  0.00139209  1.58489e9  2420.0   2900.0   0.859079  3.94439e12      3.94903e10       3.94903e9  0.149114
   2 │ 100000.0        2.15823e8  1405.74  0.00139197  1.58489e9  2420.0   2900.0   0.859081  3.9444e12       3.94903e10       3.94903e9  0.149115
   3 │      6.31887e5  2.1584e8   1405.74  0.00139134  1.5849e9   2420.0   2900.01  0.859088  3.94442e12      3.94905e10       3.94905e9  0.149122
  ⋮  │       ⋮             ⋮         ⋮         ⋮           ⋮         ⋮        ⋮        ⋮          ⋮             ⋮               ⋮            ⋮
 456 │      2.98793e9  2.21924e8  1401.0   0.00139538  1.58945e9  2421.59  2901.91  0.856613  3.97744e12      3.98211e10       3.98211e9  0.17412
 457 │      2.99793e9  2.22224e8  1401.01  0.00138474  1.58951e9  2421.66  2902.0   0.856745  3.97782e12      3.9825e10        3.9825e9   0.174241
 458 │      3.0e9      2.22286e8  1401.01  0.00138255  1.58952e9  2421.68  2902.01  0.856772  3.9779e12       3.98258e10       3.98258e9  0.174266
                                                                                                                                   452 rows omitted
```
The output directory specified by `output_dirname` contains the generated files.

### Example 3: Performing Multi-Dataset Computations
This function also allows you to perform computations for all combinations of input parameters(`log_volume_km3`, `InitialConc_H2O`, `InitialConc_CO2`, `log_vfr`, `depth`) specified by multiple datasets provided as vectors.
```julia
julia> composition = Mafic()
julia> end_time = 3e9
julia> log_volume_km3 = 0.2
julia> InitialConc_H2O = [0.01, 0.02]
julia> InitialConc_CO2 = 0.001
julia> log_vfr = -3.3
julia> depth = [7e3, 8e3]
julia> output_dirname = "MyDirname"

julia> chamber(composition, end_time, log_volume_km3, InitialConc_H2O, InitialConc_CO2, log_vfr, depth, output_dirname)
Output path: YOUR_PATH\MyDirname\vol0.2_h2o0.01_gas0.001_vfr-3.3_dep7000.0
Output path: YOUR_PATH\MyDirname\vol0.2_h2o0.01_gas0.001_vfr-3.3_dep8000.0
Output path: YOUR_PATH\MyDirname\vol0.2_h2o0.02_gas0.001_vfr-3.3_dep7000.0
Output path: YOUR_PATH\MyDirname\vol0.2_h2o0.02_gas0.001_vfr-3.3_dep8000.0
"MyDirname"
```

## Notebook for Google Colaboratory
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/brown-ccv/Chamber.jl/blob/master/notebooks/notebook_for_chamber.ipynb)

## API Documentation

API documentation for Chamber.jl can be found [here][docs-stable-url].

## References
- W. Degruyter and C. Huber. A model for eruption frequency of upper crustal silicic magma chambers. Earth Planet. Sci. Lett. (2014).
- DifferentialEquations.jl. Available at https://github.com/SciML/DifferentialEquations.jl.

[docs-dev-url]: https://brown-ccv.github.io/Chamber.jl/dev
[docs-stable-url]: https://brown-ccv.github.io/Chamber.jl/stable/
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg

