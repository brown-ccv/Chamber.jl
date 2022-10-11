function make_param()
    param = Dict{Any,Any}([])
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

    param["Q_out_old"] = 0

    param["dP_lit_dt"] = 0
    param["dP_lit_dt_0"] = param["dP_lit_dt"]
    param["P_lit_drop_max"] = 9e6

    # set the mass inflow rate
    param["Mdot_out_pass"] = 10000
    return param
end

function make_param_saved_var()
    return Dict{Any,Any}(
        "maxTime" => 0,
        "lengthTime" => 0,
        "switch_Tprofile" => 0)
end

function make_sw()
        # eruption/cooling_module/viscous_relaxation control
        return Dict{Any, Any}(
            "heat_cond" => 1,   # switch cooling module on/off
            "visc_relax" => 1,   # switch viscous relaxation on/off
            "eruption" => 0
    )
end

function make_param_IC_Finder()
    # IC Finder parameters
    param_IC_Finder = Dict{Any, Any}([])
    param_IC_Finder["max_count"] = 100
    param_IC_Finder["min_eps_g"] = 1e-10
    param_IC_Finder["eps_g_guess_ini"] = 1e-2
    param_IC_Finder["X_co2_guess_ini"] = 0.2
    param_IC_Finder["fraction"] = 0.2
    param_IC_Finder["delta_X_co2"] = 1e-2
    return param_IC_Finder
end