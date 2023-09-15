module Chamber
using CSV
using DataFrames
using Dates
using DifferentialEquations
using LinearAlgebra
using Memoize
using Parameters
using Plots
using Roots
using SpecialFunctions
using Sundials
using TimerOutputs
using LaTeXStrings
using Logging
include("utils.jl")
include("./glq_points_weights.jl")
include("./initial-utils.jl")
include("./initial.jl")
include("./heat.jl")
include("./write_csv.jl")
include("./plot_figs.jl")
include("./boundary_condition.jl")
include("./runcode-func.jl")
include("./ic_finder-utils.jl")
include("./ic_finder.jl")
include("./utils-matrix.jl")
include("../scripts/runcode.jl")
export @memoize,
    Silicic,
    Mafic,
    QNDF,
    FBDF,
    Rodas4,
    Tsit5,
    KenCarp4,
    CVODE_BDF,
    Rosenbrock23,
    TimerOutput,
    @timeit,
    reset_timer!,
    print_timer,
    print_timer_log,
    SimpleLogger,
    global_logger,
    with_logger,
    get_timer,
    SciMLBase,
    VectorContinuousCallback,
    ODEProblem,
    solve,
    get_timestamp,
    OdeSetting,
    rheol_composition_dict,
    rheol_dict,
    DataFrame,
    eos_g,
    eos_g_rho_g,
    dX_dxdydz,
    parameters_melting_curve,
    crystal_fraction,
    crystal_fraction_eps_x,
    exsolve,
    exsolve_meq,
    find_liq,
    gas_heat_capacity,
    Co2PartitionCoeff,
    boundary_conditions_new,
    heat_conduction_chamber_profileCH,
    water,
    dwater_dx,
    solve_NR,
    exsolve3,
    mco2_dissolved_sat,
    meq_water,
    get_phase,
    solve_X_co2,
    get_eps_g,
    IC_Finder,
    GLQ_points_weights_hard,
    compute_dXdP_dXdT,
    odeChamber,
    stopChamber_MT,
    affect!,
    ConstantValues,
    Param,
    ParamSaved,
    ParamICFinder,
    SW,
    EruptSaved,
    ChamberOutput,
    rho_f,
    drho_dX_f,
    rc_f,
    drc_dX_f,
    build_rho_rc,
    rho_0_f,
    build_mdot_in,
    meq_silicic,
    dmeqdT_silicic,
    dmeqdP_silicic,
    dmeqdXco2_silicic,
    meq_mafic,
    dmeqdT_mafic,
    dmeqdP_mafic,
    dmeqdXco2_mafic,
    C_co2_f,
    dC_co2dT_f,
    dC_co2dP_f,
    dC_co2dXco2_f,
    build_meq,
    build_co2,
    a1x_f,
    a13_f,
    a21_f,
    a22_f,
    a23_f,
    a24_f,
    a31_f,
    a32_f,
    a33_f,
    a34_f,
    a41_f,
    a42_f,
    a43_f,
    a44_f,
    b1_f,
    b2_f,
    b3_f,
    b4_f,
    build_matrix,
    rho_f,
    drho_dX_f,
    rc_f,
    drc_dX_f,
    build_rho_rc,
    write_csv,
    plot_figs,
    plot_sol_year_unit,
    latexstringtitle,
    plot_combined_fig,
    plot_dual_axis,
    plot_ϵx,
    write_ϵx_csv,
    check_for_duplicates,
    chamber

# Write your package code here.
end
