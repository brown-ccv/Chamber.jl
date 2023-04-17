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
composition = Silicic()
end_time = 3e9
log_volume_km3 = 0.2
InitialConc_H2O = 0.04
InitialConc_CO2 = 0.001
log_vfr = -3.3
depth = 8e3

# Run simulation
df = chamber(composition, end_time, log_volume_km3, InitialConc_H2O, InitialConc_CO2, log_vfr, depth)

# The generated CSV file can be found at 'out.csv'
# Figures for P+dP, T, eps_g, V, X_CO2, total_mass are also generated and can be viewed in the output directory.
```

## API Documentation
#=
NOTE: unfinished
=#

API documentation for the chamber package can be found [here](https://your-package-docs-url-here/).

## References
- W. Degruyter and C. Huber. A model for eruption frequency of upper crustal silicic magma chambers. Earth Planet. Sci. Lett. (2014).
- DifferentialEquations.jl. Available at https://github.com/SciML/DifferentialEquations.jl.
