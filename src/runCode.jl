using LinearAlgebra
using DifferentialEquations, Sundials #CVODE_BDF
using DataFrames, CSV, TimerOutputs, Dates
using Plots
# include("func.jl")
include("runCode-func.jl")
using .myfunc: eos_g, crystal_fraction_silicic, crystal_fraction_mafic, exsolve_silicic, exsolve_mafic, find_liq_silicic, find_liq_mafic, gas_heat_capacity, IC_Finder_silicic, IC_Finder_mafic, boundary_conditions_new, heat_conduction_chamber_profileCH, exsolve3_silicic, exsolve3_mafic, GLQ_points_weights_hard

"""
    chamber(composition::String, end_time::Int64, log_volume_km3::Number, range_water::Float64, range_co2::Float64, log_vfr::Float64, range_depth::Number)

The volcano eruption simulation

# Arguments

- `composition`: "silicic" or "mafic"
- `end_time`: simulation period. ex. 3e9
- `log_volume_km3`: Estimate volume of volcano chamber in log scale.
- `range_water`: water content
- `range_co2`: CO2 content
- `log_vfr`: 
- `range_depth`: Estimate depth of volcano chamber
"""
function chamber(composition::String, end_time::Number, log_volume_km3::Number, range_water::Float64, range_co2::Float64, log_vfr::Float64, range_depth::Number) # ("silicic", 1e9, 0.2, 0.04, 0.001, -3.3, 8e3)
    if !(composition in ["silicic", "mafic"])
        @error("composition should be \"silicic\" or \"mafic\", not \"$composition\"")
        return "Stop"
    end
    plott = true
    plot_matlab = false

    method = "CVODE_BDF"
    datetime = Dates.format(now(), "YYmmddHHMM")
    path = joinpath(pwd(), "$(datetime)_$composition")
    mkdir(path)
    io = open("$path/$datetime.log", "a")

    to = get_timer("share")
    @timeit to "runCode" begin
    # time
    begin_time     = 0
    # error tolerances used in ode method
    reltol        = 1e-8
    abstol        = 1e-8
    first_step    = 1e5 # 1e6 set lower
    max_step      = 1e7 # 1e8 set higher # In silicic, T will error when max_step >= 1e8. Use 1e7.
    write(io, "rtol: $reltol, atol: $abstol\n")
    write(io, "first_step: $first_step, max_step: $max_step \n")

    methods = Dict(
        "QNDF"=>QNDF(autodiff=false),
        "FBDF"=>FBDF(autodiff=false),
        "Rodas4"=> Rodas4(autodiff=false),
        "Tsit5"=>Tsit5(),
        "KenCarp4"=>KenCarp4(autodiff=false),
        "CVODE_BDF"=>CVODE_BDF(),
        "Rosenbrock23"=>Rosenbrock23(autodiff=false),
    )

    param = Dict{Any,Any}([])
    param["composition"] = composition
    param["rheol"] = "old"
    param["fluxing"] = false
    param["single_eruption"] = false

    # some constants
    param["beta_m"]  = 1e10   # bulk modulus melt (Pa)
    param["beta_x"]  = 1e10   # bulk moduluis crystals (Pa)
    param["alpha_m"] = 1e-5   # thermal expansion melt (K^-1)
    param["alpha_x"] = 1e-5   # thermal expansion crystals (K^-1)
    param["L_e"]     = 610e3  # latent heat of exsolution (J kg^-1) value used by Caricchi and Blundy (2015)

    param["mm_co2"] = 44.01e-3   # molar mass of CO2 (kg/mol)
    param["mm_h2o"] = 18.02e-3   # molar mass of H2O (kg/mol)

    # crustal properties
    param["alpha_r"] = 1e-5    # thermal expansion country rock (K^-1)
    param["beta_r"]  = 1e10    # bulk modulus country rock (Pa)
    param["kappa"]   = 1e-6    # thermal conductivity (W m^-1 K^-1)
    param["rho_r"]   = 2750    # density crust (kg/m3)
    param["c_r"]     = 1200    # heat capacity crust (J kg^-1 K-1)

    # Rheology of the crust:
    param["maxn"] = 10000
    param["GLQ_n"] = 64

    if param["composition"] == "silicic"
        rho_m0       = 2400      # initial melt density (kg/m^3)
        rho_x0       = 2600      # initial crystal density (kg/m^3)
        param["c_m"] = 1200      # specific heat of melt (J kg^-1 K-1)
        param["c_x"] = 1200      # specific heat of crystals (J kg^-1 K-1)
        param["L_m"] = 290e3     # latent heat of melting (J kg^-1) value used by Caricchi and Blundy (2015)
    elseif param["composition"] == "mafic"
        rho_m0       = 2420      # initial melt density (kg/m^3)
        rho_x0       = 2900      # initial crystal density (kg/m^3)
        param["c_m"] = 1142      # specific heat of melt (J kg^-1 K-1)
        param["c_x"] = 1160      # specific heat of crystals (J kg^-1 K-1)
        param["L_m"] = 470e3     # latent heat of melting (J kg^-1) value used by Caricchi and Blundy (2015)
    end

    if param["rheol"] == "new"
        param["A"] = 4.25e7      # material-dependent constant for viscosity law (Pa s)
        param["B"] = 8.31        # molar gas constant (J/mol/K)
        param["G"] = 141e3       # activation energy for creep (J/mol)
    elseif param["rheol"] == "old"
        param["nn"] = 1.9        # powerlaw exponent
        param["AA"] = 2e-4*(1e6)^-param["nn"] #material dependent constant (Pa^-n s^-1)
        param["G"]  = 141e3      # activation energy for creep (J/mol)
        param["M"]  = 8.31       # gas constant
    end

    param_saved_var = Dict{Any, Any}(
        "maxTime" => 0,
        "lengthTime" => 0,
        "switch_Tprofile" => 0)
    param_saved_var["storeSumk"] = zeros(param["maxn"])
    param_saved_var["storeSumk_2"] = zeros(param["maxn"])
    param_saved_var["storeSumk_old"] = zeros(param["maxn"])
    param_saved_var["storeSumk_2_old"]= zeros(param["maxn"])
    param["Q_out_old"] = 0

    range_volume_km3     = 10^log_volume_km3   # range of volume in km3
    range_radius         = 1000*(range_volume_km3/(4*pi/3))^(1/3)   # range of radius in m
    V_0  = 4*pi/3*range_radius^3   # initial volume of the chamber (m^3)
    range_vfr        = 10^log_vfr   # volume flow rate (km3/yr)  
    range_mfr        = if param["fluxing"] 0 else rho_m0*range_vfr*1e9/(3600*24*365) end

    println("composition: $composition,  method: $method")
    println("log_volume_km3: ", log_volume_km3, " range_water: ", range_water, " range_co2: ", range_co2)
    println("log_vfr: ", log_vfr, " range_depth: ", range_depth)
    write(io, "composition: $composition \n")
    write(io, "log_volume_km3: $log_volume_km3, range_water: $range_water, range_co2: $range_co2 \n")
    write(io, "log_vfr: $log_vfr, range_depth: $range_depth \n")
    write(io, "method: $method \n")

    ## initial conditions
    # depth reservoir
    depth = range_depth
    InitialConc_H2O = range_water
    InitialConc_CO2 = range_co2

    # thermal gradient
    T_surface     = 0+273   # surface temperature (K)
    T_gradient    = 32/1e3   # thermal gradient (K/m)
    Tb            = T_surface+T_gradient*depth   # background temperature crust (K)

    # lithostatic pressure
    grav_acc      = 9.81   # gravitational acceleration (m/s2)
    DP_crit       = 20e6   # critical overpressure (Pa)
    P_0           = param["rho_r"]*grav_acc*depth   # initial chamber pressure (Pa)
    P_lit         = P_0
    param["Tb"]   = Tb
    param["P_lit_0"] = P_0
    param["dP_lit_dt"] = 0
    param["dP_lit_dt_0"] = param["dP_lit_dt"]

    param["P_lit_drop_max"] = 9e6
    P0plusDP_0  = param["P_lit_0"]

    if param["single_eruption"]
        P_0 = P_0 + DP_crit
        P0plusDP_0 = P_0
    end

    ini_eps_x  = 0.15
    if param["composition"] == "silicic"
        @timeit to "find_liq_silicic" T_0 = find_liq_silicic(InitialConc_H2O, InitialConc_CO2, P_0, ini_eps_x)
    elseif param["composition"] == "mafic"
        @timeit to "find_liq_mafic" T_0 = find_liq_mafic(InitialConc_H2O, InitialConc_CO2, P_0, ini_eps_x)
    end
    T_R  = T_0+50   # Temperature of recharging magma (K)
    T_in = T_R   # Temperature of inflowing magma (K)

    if ~param["fluxing"]
        mdot_in = range_mfr
    else
        log_vfr    = -4.3   # log volume flow rate (km3/yr)
        range_vfr        = 10^log_vfr   # volume flow rate (km3/yr)
        rho_g_in         = eos_g(P_0, T_in)["rho_g"]
        range_mfr        = rho_g_in*range_vfr*1e9/(3600*24*365)
        mdot_in          = range_mfr   # place holder
        XCO2_in          = 0.8
        param["XCO2_in"] = XCO2_in
    end

    # set the mass inflow rate
    param["Mdot_in_pass"] = mdot_in
    param["Mdot_out_pass"] = 10000

    rho_g0 = eos_g(P_0, T_0)["rho_g"]   # initial gas density
    if param["composition"] == "silicic"
        @timeit to "crystal_fraction_silicic" eps_x0 = crystal_fraction_silicic(T_0, P_0, InitialConc_H2O, InitialConc_CO2)[1]
    elseif param["composition"] == "mafic"
        @timeit to "crystal_fraction_mafic" eps_x0 = crystal_fraction_mafic(T_0, P_0, InitialConc_H2O, InitialConc_CO2)[1]
    end

    eps_m0 = 1-eps_x0
    rho = rho_m0*eps_m0 + rho_x0*eps_x0
    M_tot =  V_0*rho   # Total mass, initial
    M_co2_0 = InitialConc_CO2*V_0*rho   # Total mass of CO2, initial
    M_h2o_0 = InitialConc_H2O*V_0*rho   # Total mass of H2O, initial

    # IC Finder parameters
    param_IC_Finder = Dict{Any, Any}([])
    param_IC_Finder["max_count"] = 100
    param_IC_Finder["Tol"] = if composition == "silicic" 1e-9 else 1e-8 end
    param_IC_Finder["min_eps_g"] = 1e-10
    param_IC_Finder["eps_g_guess_ini"] = 1e-2
    param_IC_Finder["X_co2_guess_ini"] = 0.2
    param_IC_Finder["fraction"] = 0.2
    param_IC_Finder["delta_X_co2"] = 1e-2
    write(io, "IC_Finder parameters: $(param_IC_Finder)\n")

    if param["composition"] == "silicic"
        @timeit to "IC_Finder_silicic" eps_g0, X_co20, C_co2, phase = IC_Finder_silicic(M_h2o_0, M_co2_0, M_tot, P_0, T_0, V_0, rho_m0, param["mm_co2"], param["mm_h2o"], param_IC_Finder)
    elseif param["composition"] == "mafic"
        @timeit to "IC_Finder_mafic" eps_g0, X_co20, C_co2, phase = IC_Finder_mafic(M_h2o_0, M_co2_0, M_tot, P_0, T_0, V_0, rho_m0, rho_x0, param["mm_co2"], param["mm_h2o"], param_IC_Finder)
    end
    println("IC_Finder done: [eps_g0, X_co20, C_co2] = [$eps_g0, $X_co20, $C_co2]")
    println("phase: ", phase)
    write(io, "IC_Finder done: $eps_g0, $X_co20, $C_co2, $phase\n")
    param_saved_var["phase"] = phase

    # update initial crystal volume fraction
    if param["composition"] == "silicic"
        @timeit to "crystal_fraction_silicic" eps_x0 = crystal_fraction_silicic(T_0, P_0, InitialConc_H2O, InitialConc_CO2)[1]
    elseif param["composition"] == "mafic"
        @timeit to "crystal_fraction_mafic" eps_x0 = crystal_fraction_mafic(T_0, P_0, InitialConc_H2O, InitialConc_CO2)[1]
    end

    # update initial melt volume fraction
    eps_m0 = 1-eps_x0-eps_g0
    # update initial bulk density (kg/m^3)
    rho_0  = (1-eps_g0-eps_x0)*rho_m0 + eps_g0*rho_g0 + eps_x0*rho_x0

    # update solubility
    if phase == 2
        m_eq0 = M_h2o_0/((1-eps_x0)*rho_m0*V_0)
        C_co20 = M_co2_0/((1-eps_x0)*rho_m0*V_0)
        X_co20 = 0
    else
        if param["composition"] == "silicic"
            @timeit to "exsolve_silicic" ans_es = exsolve_silicic(P_0, T_0, X_co20)
            m_eq0, C_co20 = ans_es[1], ans_es[5]
        elseif param["composition"] == "mafic"
            @timeit to "exsolve_mafic" ans_es = exsolve_mafic(P_0, T_0, X_co20)
            m_eq0, C_co20 = ans_es[1], ans_es[5]
        end
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

    sw = Dict{Any, Any}(
        "heat_cond" => 1,   # switch cooling module on/off
        "visc_relax" => 1,   # switch viscous relaxation on/off
        "eruption" => 0
    )
    if param["single_eruption"]
        sw["eruption"] = 1
    end

    # initialize vector to store quantities
    storeTime = [0]
    storeTemp = [T_0]

    # Initial temperature and viscosity profile around chamber
    @timeit to "GLQ_points_weights_hard" quadpts, weights = GLQ_points_weights_hard(param["GLQ_n"])
    cc    = 10*range_radius   #outer radius for heat conduction (m)
    b     = range_radius + cc
    quadpts_r = (b-range_radius)/2*quadpts .+ (range_radius+b)/2
    param_saved_var["storeTime"] = storeTime
    param_saved_var["storeTemp"] = storeTemp

    @timeit to "heat_conduction_chamber_profileCH" Trt = heat_conduction_chamber_profileCH(param["maxn"], range_radius, cc, quadpts_r, param["kappa"], param["Tb"], param_saved_var)
    if param["rheol"] == "new"
        A = param["A"]
        B = param["B"]
        G = param["G"]
        eta_rt     = A*exp(G/B/Trt)
    elseif param["rheol"] == "old"
        nn = param["nn"]
        AA = param["AA"]
        G = param["G"]
        M = param["M"]
        dev_stress = DP_crit
        eta_rt     = (dev_stress^(1-nn)/AA)*exp(G/M/Trt)
    end

    param["P_lit"] = P_lit
    param["DP_crit"] = DP_crit

    tspan    = (begin_time, end_time)
    IC       = [P0plusDP_0, T_0, eps_g0, V_0, rho_m0, rho_x0, X_co20, tot_Mass_0, tot_Mass_H2O_0, tot_Mass_CO2_0]
    println("tspan: ", tspan)
    println("IC: ", IC)
    write(io, "IC: $IC\n")
    write(io, "  t, P+dP, eps_g, X_co2 : \n")
    cb = VectorContinuousCallback(stopChamber_MT, affect!, 8, rootfind=SciMLBase.RightRootFind)
    prob = ODEProblem(odeChamber,IC,tspan,param)
    @timeit to "solve_ODE" sol = solve(prob, methods[method], callback=cb, reltol=reltol, abstol=abstol, dt=first_step, dtmax=max_step)

    end
    println(to)
    write(io, "$to\n")

    df = DataFrame(sol)
    number_of_data = length(sol)
    println("number_of_data: $number_of_data")
    write(io, "number_of_data: $number_of_data\n")
    CSV.write("$path/out.csv", df)
    close(io)

    ENV["GKSwstype"] = "100"  # magic environmental variable for Plots
    if plott == true && plot_matlab === true
        ## plot matlab result
        matlab = joinpath(pwd(),"/home/calliehsu/Python_code/Mafic/matlab_mafic_wa_1e-2_co_8e-4_step_1e5_1e7.csv")  #matlab_mafic_wa_1e-2_co_1e-3 #matlab_mafic #matlab_silicic  # matlab_silicic_sat_005_00001_initialsteps_1e3_strange # matlab_silicic_saturated
        df2 = CSV.File(matlab) |> DataFrame
        tt = df2[!,1][df2[!,1].<=end_time]
        p1 = plot(tt, df2[!,2][1:length(tt)], xaxis="Time (t)", yaxis="P+dP", label="Matlab", linewidth=2, marker=(:circle,3,Plots.stroke(1,:gray)))
        p2 = plot(tt, df2[!,3][1:length(tt)], xaxis="Time (t)",yaxis="T", label="Matlab", linewidth=2, marker=(:circle,3,Plots.stroke(1,:gray)))
        p3 = plot(tt, df2[!,4][1:length(tt)], xaxis="Time (t)",yaxis="eps_g", label="Matlab", linewidth=2, marker=(:circle,3,Plots.stroke(1,:gray)))
        p4 = plot(tt, df2[!,5][1:length(tt)], xaxis="Time (t)",yaxis="V", label="Matlab", linewidth=2, marker=(:circle,3,Plots.stroke(1,:gray)))
        p5 = plot(tt, df2[!,8][1:length(tt)], xaxis="Time (t)",yaxis="X_CO2", label="Matlab", linewidth=2, marker=(:circle,3,Plots.stroke(1,:gray)))
        p6 = plot(tt, df2[!,9][1:length(tt)], xaxis="Time (t)",yaxis="tot_Mass", label="Matlab", linewidth=2, marker=(:circle,3,Plots.stroke(1,:gray)))
        plot!(p1, df[!,1], df[!,2], label="Julia $method", linewidth=2, marker=(:x,3,Plots.stroke(2)))
        plot!(p2, df[!,1], df[!,3], label="Julia $method", linewidth=2, marker=(:x,3,Plots.stroke(2)))
        plot!(p3, df[!,1], df[!,4], label="Julia $method", linewidth=2, marker=(:x,3,Plots.stroke(2)))
        plot!(p4, df[!,1], df[!,5], label="Julia $method", linewidth=2, marker=(:x,3,Plots.stroke(2)))
        plot!(p5, df[!,1], df[!,8], label="Julia $method", linewidth=2, marker=(:x,3,Plots.stroke(2)))
        plot!(p6, df[!,1], df[!,9], label="Julia $method", linewidth=2, marker=(:x,3,Plots.stroke(2)))
    elseif plott == true
        p1 = plot(df[!,1], df[!,2], xaxis="Time (t)",yaxis="P+dP", label="Julia $method", linewidth=2, marker=(:x,3,Plots.stroke(2)))
        p2 = plot(df[!,1], df[!,3], xaxis="Time (t)",yaxis="T", label="Julia $method", linewidth=2, marker=(:x,3,Plots.stroke(2)))
        p3 = plot(df[!,1], df[!,4], xaxis="Time (t)",yaxis="eps_g", label="Julia $method", linewidth=2, marker=(:x,3,Plots.stroke(2)))
        p4 = plot(df[!,1], df[!,5], xaxis="Time (t)",yaxis="V", label="Julia $method", linewidth=2, marker=(:x,3,Plots.stroke(2)))
        p5 = plot(df[!,1], df[!,8], xaxis="Time (t)",yaxis="X_CO2", label="Julia $method", linewidth=2, marker=(:x,3,Plots.stroke(2)))
        p6 = plot(df[!,1], df[!,9], xaxis="Time (t)",yaxis="tot_Mass", label="Julia $method", linewidth=2, marker=(:x,3,Plots.stroke(2)))
    end

    if plott == true
        savefig(p1, "$path/1_P.png")
        savefig(p2, "$path/2_T.png")
        savefig(p3, "$path/3_eps_g.png")
        savefig(p4, "$path/4_V.png")
        savefig(p5, "$path/5_X_CO2.png")
        savefig(p6, "$path/6_tot_Mass.png")
    end

    println(".. Done!")
end
