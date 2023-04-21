[![Build Status](https://github.com/CallieHsu/Chamber.jl/actions/workflows/test.yml/badge.svg?event=push)](https://github.com/CallieHsu/Chamber.jl/actions/workflows/test.yml??query=branch%3Amain)

# Chamber.jl
`Chamber` is a Julia package for simulating the eruption of a volcano using a model for the frequency of eruptions of upper crustal magma chambers based on Degruyter and Huber (2014). The package generates a CSV file and figures with the solution data of various variables over time.

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
Here's an example of how to use the chamber function:

```julia
using Chamber

# Define simulation parameters
composition = Silicic() # chose here between Mafic() – basalt composition (rift) – rhyolite composition (arc magma)
end_time = 3e9 # maximum magma chamber evolution duration in seconds
log_volume_km3 = 0.2 # initial chamber volume is 10^(log_volume_km3) in km3
InitialConc_H2O = 0.04 # initial weight fraction of water in the magma (exsolved+dissolved)
InitialConc_CO2 = 0.001 # initial weight fraction of CO2 in the magma (exsolved+dissolved)
log_vfr = -3.3 # magma recharge rate [km3/yr] = 10^(log_vfr)
depth = 8e3 # depth of the magma chamber in meters

# Run simulation
df = chamber(composition, end_time, log_volume_km3, InitialConc_H2O, InitialConc_CO2, log_vfr, depth)

# Two generated CSV files 
# File 1: 'out.csv' in folder X, it is a comma-delimited ascii file organized as follows : 1st col time (sec), 2nd col P+dP (Pa), 3rd col T (K), 4th col eps_g, 5th col V (m³), 6th col rho_m (kg/m³), 7th col rho_x (kg/m³), 8th col X_CO2, 9th col total_mass (kg), 10th col total_mass_H2O (kg), 11th col total_mass_CO2 (kg)
# File 2: 'eruptions.csv'  in folder X, also a comma-delimited ascii file. It is organized as follows: 1st col time of eruption (sec), 2nd col. Duration of eruption (sec), 3rd col. Mass erupted (kg), 4th col. Volume erupted (km3)
# Figures for P+dP(t), T(t), eps_g(t), V(t), X_CO2(t), total_mass(t) are generated automatically and can be viewed in the output directory (Name provided by new variable).

```

## API Documentation
#=
NOTE: unfinished
=#

API documentation for the chamber package can be found [here](https://your-package-docs-url-here/).

## References
- W. Degruyter and C. Huber. A model for eruption frequency of upper crustal silicic magma chambers. Earth Planet. Sci. Lett. (2014).
- DifferentialEquations.jl. Available at https://github.com/SciML/DifferentialEquations.jl.
