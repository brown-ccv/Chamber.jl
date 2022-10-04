module Chamber
using CSV
using DataFrames
using Dates
using DifferentialEquations
using LinearAlgebra
using Plots
using Roots
using SpecialFunctions
using Sundials
using TimerOutputs
include("./GLQ_points_weights.jl")
include("./initial.jl")
include("./heat.jl")
include("./write_csv.jl")
include("./plot_figs.jl")
include("./boundary_condition.jl")
include("./runCode-func.jl")
include("./IC_finder.jl")
export eos_g,
       crystal_fraction_mafic,
       exsolve_silicic,
       exsolve_mafic,
       find_liq_silicic,
       find_liq_mafic,
       gas_heat_capacity,
       IC_Finder_silicic,
       IC_Finder_mafic,
       boundary_conditions_new,
       heat_conduction_chamber_profileCH,
       exsolve3_silicic,
       exsolve3_mafic,
       GLQ_points_weights_hard

# Write your package code here.

end
