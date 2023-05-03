[![Build Status](https://github.com/brown-ccv/Chamber.jl/actions/workflows/test.yml/badge.svg)](https://github.com/brown-ccv/Chamber.jl/actions?query=workflows/test)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]

# Chamber.jl
`Chamber.jl` is a Julia package for simulating the eruption of a volcano using a model for the frequency of eruptions of upper crustal magma chambers based on Degruyter and Huber (2014). The package generates a CSV file and figures with the solution data for various variables over time.

## Installation
#=
NOTE: unfinished
=#

To install the `Chamber` package, open Julia and use the package manager:

```julia
using Pkg
Pkg.add("Chamber")
```

## Usage
The main function in the `Chamber` package is `chamber`, which simulates the evolution of a magma chamber over time and returns a `DataFrame` with the solution data. The function takes the following arguments:
```julia
chamber(composition, end_time, log_volume_km3, InitialConc_H2O, InitialConc_CO2, log_vfr, depth, output_dirname)
```
### Arguments
- `composition`: The magma composition. Use `Silicic()` for rhyolite composition (arc magma) or `Mafic()` for basalt composition (rift).
- `end_time`: Maximum magma chamber evolution duration in seconds.
- `log_volume_km3`: The initial volume of the chamber in logarithmic scale. The actual initial chamber volume is calculated as 10^(log_volume_km3) in km³.
- `InitialConc_H2O`: The initial weight fraction of water in the magma (exsolved + dissolved).
- `InitialConc_H2O`: The initial weight fraction of CO₂ in the magma (exsolved + dissolved).
- `log_vfr`: Magma recharge rate in km³/yr calculated as 10^(`log_vfr`).
- `depth`: Depth of the magma chamber in meters.
- `output_dirname`: (Optional) Name of the output directory. Defaults to current timestamp.

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
- `total_mass`: Total mess of magma chamber (kg).
- `total_mass_H2O`: Total mess of water in the magma (kg).
- `total_mass_CO2`: Total mass of CO₂ in the magma (kg).

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
470×11 DataFrame
 Row │ time            P+dP       T        eps_g       V          rho_m    rho_x    X_CO2     total_mass  total_mass_H2O  total_mass_CO2
     │ Float64         Float64    Float64  Float64     Float64    Float64  Float64  Float64   Float64     Float64         Float64
─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │      0.0        2.1582e8   1046.71  0.00558159  1.58489e9  2400.0   2600.0   0.346043  3.83394e12      1.54039e11       3.85098e9
   2 │ 100000.0        2.15824e8  1046.71  0.00558131  1.58489e9  2400.0   2600.0   0.34605   3.83395e12      1.54039e11       3.85099e9
   3 │      7.62058e5  2.15848e8  1046.71  0.0055795   1.5849e9   2400.01  2600.01  0.346101  3.83397e12      1.5404e11        3.85101e9
   4 │      2.83015e6  2.15924e8  1046.71  0.00557383  1.58491e9  2400.03  2600.03  0.346259  3.83405e12      1.54044e11       3.85109e9
   5 │      4.89825e6  2.16e8     1046.71  0.00556818  1.58492e9  2400.04  2600.05  0.346418  3.83413e12      1.54047e11       3.85117e9
  ⋮  │       ⋮             ⋮         ⋮         ⋮           ⋮         ⋮        ⋮        ⋮          ⋮             ⋮               ⋮
 467 │      2.9764e9   2.23039e8  1043.45  0.00561742  1.58776e9  2401.81  2601.96  0.337852  3.8504e12       1.547e11         3.86751e9
 468 │      2.9864e9   2.23399e8  1043.44  0.00559152  1.58782e9  2401.9   2602.06  0.338571  3.85078e12      1.54716e11       3.86789e9
 469 │      2.9964e9   2.2376e8   1043.43  0.00556579  1.58788e9  2401.99  2602.15  0.339287  3.85116e12      1.54731e11       3.86828e9
 470 │      3.0e9      2.23889e8  1043.43  0.00555656  1.5879e9   2402.02  2602.18  0.339545  3.8513e12       1.54737e11       3.86841e9
                                                                                                                         461 rows omitted
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
458×11 DataFrame
 Row │ time            P+dP       T        eps_g       V          rho_m    rho_x    X_CO2     total_mass  total_mass_H2O  total_mass_CO2 
     │ Float64         Float64    Float64  Float64     Float64    Float64  Float64  Float64   Float64     Float64         Float64        
─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │      0.0        2.1582e8   1405.74  0.00139209  1.58489e9  2420.0   2900.0   0.859079  3.94439e12      3.94903e10       3.94903e9
   2 │ 100000.0        2.15823e8  1405.74  0.00139197  1.58489e9  2420.0   2900.0   0.859081  3.9444e12       3.94903e10       3.94903e9
  ⋮  │       ⋮             ⋮         ⋮         ⋮           ⋮         ⋮        ⋮        ⋮          ⋮             ⋮               ⋮
 457 │      2.99798e9  2.22224e8  1401.01  0.00138553  1.58951e9  2421.66  2901.99  0.856719  3.97782e12      3.98249e10       3.98249e9
 458 │      3.0e9      2.22284e8  1401.01  0.00138339  1.58952e9  2421.68  2902.01  0.856746  3.9779e12       3.98257e10       3.98257e9
                                                                                                                         454 rows omitted
```
The output directory specified by `output_dirname` contains the generated files.

## API Documentation

API documentation for Chamber.jl can be found [here][docs-dev-url].

## References
- W. Degruyter and C. Huber. A model for eruption frequency of upper crustal silicic magma chambers. Earth Planet. Sci. Lett. (2014).
- DifferentialEquations.jl. Available at https://github.com/SciML/DifferentialEquations.jl.

[docs-dev-url]: https://brown-ccv.github.io/Chamber.jl/
