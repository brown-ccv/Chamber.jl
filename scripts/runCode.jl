using Chamber
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
using Logging
include("./solver_methods.jl")
include("./utils.jl")

# thermal gradient
T_surface     = 0+273   # surface temperature (K)
T_gradient    = 32/1e3   # thermal gradient (K/m)

# lithostatic pressure
grav_acc      = 9.81   # gravitational acceleration (m/s2)
DP_crit       = 20e6   # critical overpressure (Pa)
"""
    chamber(composition::String, end_time::Int64, log_volume_km3::Number, InitialConc_H2O::Float64, InitialConc_CO2::Float64, log_vfr::Float64, depth::Number)

The volcano eruption simulation

# Arguments

- `composition`: "silicic" or "mafic"
- `end_time`: simulation period. ex. 3e9
- `log_volume_km3`: Estimate volume of volcano chamber in log scale.
- `InitialConc_H2O`: water content
- `InitialConc_CO2`: CO2 content
- `log_vfr`: 
- `depth`: Estimate depth of volcano chamber
"""
function chamber(composition::String, end_time::Number, log_volume_km3::Number, InitialConc_H2O::Float64, InitialConc_CO2::Float64, 
    log_vfr::Float64, depth::Number, methods::Dict=methods, method::String="Tsit5", errtol=ErrTol(), ini_eps_x::Float64=0.15, rheol::String="old") # ("silicic", 1e9, 0.2, 0.04, 0.001, -3.3, 8e3)
    if !(composition in ["silicic", "mafic"])
        @error("composition should be \"silicic\" or \"mafic\", not \"$composition\"")
        return "Stop"
    end

    datetime = get_timestamp()
    path = joinpath(pwd(), "$(datetime)_$composition")
    mkdir(path)
    io = open("$path/$datetime.log", "w+")
    logger = SimpleLogger(io)
    global_logger(logger)

    param = make_param(rheol)
    param_saved_var = make_param_saved_var()
    param_IC_Finder = make_param_IC_Finder()
    sw = make_sw()
    param["composition"] = composition
    c = composition_dict[composition]
    param["c_m"] = c.c_m
    param["c_x"] = c.c_x
    param["L_m"] = c.L_m

    param_saved_var["storeSumk"] = zeros(param["maxn"])
    param_saved_var["storeSumk_2"] = zeros(param["maxn"])
    param_saved_var["storeSumk_old"] = zeros(param["maxn"])
    param_saved_var["storeSumk_2_old"]= zeros(param["maxn"])

    volume_km3    = 10^log_volume_km3                 # range of volume in km3
    range_radius  = 1000*(volume_km3/(4*pi/3))^(1/3)  # range of radius in m
    V_0           = 4*pi/3*range_radius^3             # initial volume of the chamber (m^3)

    # thermal gradient
    Tb            = T_surface+T_gradient*depth        # background temperature crust (K)

    # lithostatic pressure
    P_0           = param["rho_r"]*grav_acc*depth   # initial chamber pressure (Pa)
    P_lit         = P_0
    param["Tb"]   = Tb
    param["P_lit_0"] = P_0
    P0plusDP_0  = P_0

    if param["single_eruption"]
        P_0 = P_0 + DP_crit
        P0plusDP_0 = P_0
    end

    if composition == "silicic"
        T_0 = find_liq_silicic(InitialConc_H2O, InitialConc_CO2, P_0, ini_eps_x)
    elseif composition == "mafic"
        T_0 = find_liq_mafic(InitialConc_H2O, InitialConc_CO2, P_0, ini_eps_x)
    end
    T_R  = T_0 + 50   # Temperature of recharging magma (K)
    T_in = T_R        # Temperature of inflowing magma (K)
    param["T_in"] = T_in

    if ~param["fluxing"]  # fluxing = false
        range_vfr        = 10^log_vfr   # volume flow rate (km3/yr)  
        range_mfr        = c.rho_m0*range_vfr*1e9/(3600*24*365)
        mdot_in = range_mfr
    else  # fluxing = true
        log_vfr          = -4.3   # log volume flow rate (km3/yr)
        range_vfr        = 10^log_vfr   # volume flow rate (km3/yr)
        rho_g_in         = eos_g(P_0, T_in)["rho_g"]
        range_mfr        = rho_g_in*range_vfr*1e9/(3600*24*365)
        mdot_in          = range_mfr   # place holder
        param["XCO2_in"] = 0.8
    end

    # set the mass inflow rate
    param["Mdot_in_pass"] = mdot_in

    rho_g0 = eos_g(P_0, T_0)["rho_g"]   # initial gas density
    if composition == "silicic"
        eps_x0 = crystal_fraction_silicic(T_0, P_0, InitialConc_H2O, InitialConc_CO2)[1]
    elseif composition == "mafic"
        eps_x0 = crystal_fraction_mafic(T_0, P_0, InitialConc_H2O, InitialConc_CO2)[1]
    end

    eps_m0 = 1-eps_x0
    rho = c.rho_m0*eps_m0 + c.rho_x0*eps_x0
    M_tot =  V_0*rho   # Total mass, initial
    M_co2_0 = InitialConc_CO2*V_0*rho   # Total mass of CO2, initial
    M_h2o_0 = InitialConc_H2O*V_0*rho   # Total mass of H2O, initial

    # IC Finder parameters
    param_IC_Finder["Tol"] = if composition == "silicic" 1e-9 else 1e-8 end

    if composition == "silicic"
        eps_g0, X_co20, C_co2, phase = IC_Finder_silicic(M_h2o_0, M_co2_0, M_tot, P_0, T_0, V_0, c.rho_m0, param["mm_co2"], param["mm_h2o"], param_IC_Finder)
    elseif composition == "mafic"
        eps_g0, X_co20, C_co2, phase = IC_Finder_mafic(M_h2o_0, M_co2_0, M_tot, P_0, T_0, V_0, c.rho_m0, param["mm_co2"], param["mm_h2o"], param_IC_Finder)
    end
    println("IC_Finder done: [eps_g0, X_co20, C_co2] = [$eps_g0, $X_co20, $C_co2]")
    println("phase: ", phase)
    @info("IC_Finder done: $eps_g0, $X_co20, $C_co2, $phase")
    param_saved_var["phase"] = phase

    # update initial crystal volume fraction
    if composition == "silicic"
        eps_x0 = crystal_fraction_silicic(T_0, P_0, InitialConc_H2O, InitialConc_CO2)[1]
    elseif composition == "mafic"
        eps_x0 = crystal_fraction_mafic(T_0, P_0, InitialConc_H2O, InitialConc_CO2)[1]
    end

    # update initial bulk density (kg/m^3)
    rho_0  = (1-eps_g0-eps_x0)*c.rho_m0 + eps_g0*rho_g0 + eps_x0*c.rho_x0

    # update solubility
    if phase == 2
        X_co20 = 0
    end

    # Calculate the water content (concentration) for inflowing magma contents
    tot_h2o_frac_in = M_h2o_0/(rho_0*V_0)   # CHANGE MASS FRACTION FROM XCO2_IN
    tot_co2_frac_in = M_co2_0/(rho_0*V_0)
    param["tot_h2o_frac_in"] = tot_h2o_frac_in
    param["tot_co2_frac_in"] = tot_co2_frac_in
    if param["fluxing"]
        tot_co2_frac_in = param["XCO2_in"]*param["mm_co2"]/(param["XCO2_in"]*param["mm_co2"]+(1-param["XCO2_in"])*param["mm_h2o"])
        tot_h2o_frac_in = 1-tot_co2_frac_in
    end

    tot_Mass_0 = V_0*rho_0
    tot_Mass_H2O_0 = M_h2o_0
    tot_Mass_CO2_0 = M_co2_0

    if param["single_eruption"]
        sw["eruption"] = 1
    end

    # initialize vector to store quantities
    param_saved_var["storeTime"] = [0]
    param_saved_var["storeTemp"] = [T_0]

    param["P_lit"] = P_lit
    param["DP_crit"] = DP_crit

    println(param)
    println(param_IC_Finder)
    @info("sw: $(sw)")
    @info("IC_Finder parameters: $(param_IC_Finder)")
    @info("params: $(param)")

    tspan    = (0, end_time)
    IC       = [P0plusDP_0, T_0, eps_g0, V_0, c.rho_m0, c.rho_x0, X_co20, tot_Mass_0, tot_Mass_H2O_0, tot_Mass_CO2_0]
    println("tspan: ", tspan)
    println("IC: ", IC)
    @info("IC: $IC")
    close(io)
    cb = VectorContinuousCallback(stopChamber_MT, affect!, 8, rootfind=SciMLBase.RightRootFind)
    prob = ODEProblem(odeChamber,IC,tspan,param)
    sol = solve(prob, methods[method], callback=cb, reltol=errtol.reltol, abstol=errtol.abstol, dt=errtol.first_step, dtmax=errtol.max_step)

    write_csv(sol, path)
    plot_figs("$path/out.csv", path)

    println(".. Done!")
end

chamber("silicic", 1e9, 0.2, 0.04, 0.001, -3.3, 8e3)
