module Chamber
using CSV
using DataFrames
using Dates
using DifferentialEquations
using LinearAlgebra
using Parameters
using Plots
using Roots
using SpecialFunctions
using Sundials
using TimerOutputs
include("utils.jl")
include("./Data/parameters.jl")
include("./GLQ_points_weights.jl")
include("./initial-utils.jl")
include("./initial.jl")
include("./heat.jl")
include("./write_csv.jl")
include("./plot_figs.jl")
include("./boundary_condition.jl")
include("./runCode-func.jl")
include("./IC_finder.jl")
include("./utils-matrix.jl")
export get_timestamp, OdeSetting, rheol_composition_dict, rheol_dict, 
       eos_g,
       eos_g_rho_g,
       crystal_fraction,
       crystal_fraction_eps_x,
       exsolve,
       exsolve_meq,
       find_liq,
       gas_heat_capacity,
       IC_Finder_silicic,
       IC_Finder_mafic,
       boundary_conditions_new,
       heat_conduction_chamber_profileCH,
       exsolve3_silicic,
       exsolve3_mafic,
       GLQ_points_weights_hard,
       odeChamber,
       stopChamber_MT,
       affect!,
       make_param,
       make_param_saved_var,
       make_sw,
       make_param_IC_Finder,
       SW,
       rho_f, drho_dX_f, rc_f, drc_dX_f, build_rho_rc, rho_0_f,
       meq_silicic, dmeqdT_silicic, dmeqdP_silicic, dmeqdXco2_silicic,
       meq_mafic, dmeqdT_mafic, dmeqdP_mafic, dmeqdXco2_mafic,
       C_co2_f, dC_co2dT_f, dC_co2dP_f, dC_co2dXco2_f,
       build_meq_silicic, build_meq_mafic, build_co2,
       a1x_f, a13_f, a21_f, a22_f, a23_f, a24_f, a31_f, a32_f, a33_f, a34_f, a41_f, a42_f, a43_f, a44_f,
       b1_f, b2_f, b3_f, b4_f,
       build_matrix,
       rho_f, drho_dX_f, rc_f, drc_dX_f, build_rho_rc,
       write_csv, plot_figs

# Write your package code here.
end
