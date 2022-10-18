function make_param(composition, rheol)
    param = Dict{Any,Any}([])
    param["composition"]     = composition
    param["rheol"]           = rheol
    param["fluxing"]         = p.fluxing
    param["single_eruption"] = p.single_eruption

    # some constants
    param["beta_m"]  = p.beta_m   # bulk modulus melt (Pa)
    param["beta_x"]  = p.beta_x   # bulk moduluis crystals (Pa)
    param["alpha_m"] = p.alpha_m  # thermal expansion melt (K^-1)
    param["alpha_x"] = p.alpha_x  # thermal expansion crystals (K^-1)
    param["L_e"]     = p.L_e      # latent heat of exsolution (J kg^-1) value used by Caricchi and Blundy (2015)
    param["mm_co2"]  = p.mm_co2   # molar mass of CO2 (kg/mol)
    param["mm_h2o"]  = p.mm_h2o   # molar mass of H2O (kg/mol)

    # crustal properties
    param["alpha_r"] = p.alpha_r    # thermal expansion country rock (K^-1)
    param["beta_r"]  = p.beta_r    # bulk modulus country rock (Pa)
    param["kappa"]   = p.kappa    # thermal conductivity (W m^-1 K^-1)
    param["rho_r"]   = p.rho_r    # density crust (kg/m3)
    param["c_r"]     = p.c_r    # heat capacity crust (J kg^-1 K-1)

    # Rheology of the crust:
    param["maxn"]  = p.maxn
    param["GLQ_n"] = p.GLQ_n

    param["Q_out_old"] = p.Q_out_old
    param["dP_lit_dt"] = p.dP_lit_dt
    param["dP_lit_dt_0"] = p.dP_lit_dt_0
    param["P_lit_drop_max"] = p.P_lit_drop_max

    # set the mass inflow rate
    param["Mdot_out_pass"] = p.Mdot_out_pass

    rc = rheol_composition_dict[composition]
    param["c_m"] = rc.c_m
    param["c_x"] = rc.c_x
    param["L_m"] = rc.L_m

    r = rheol_dict[rheol]
    if rheol == "new"
        param["A"] = r.A        # material-dependent constant for viscosity law (Pa s)
        param["B"] = r.B        # molar gas constant (J/mol/K)
        param["G"] = r.G        # activation energy for creep (J/mol)
    elseif rheol == "old"
        param["nn"] = r.nn      # powerlaw exponent
        param["AA"] = r.AA      #material dependent constant (Pa^-n s^-1)
        param["G"]  = r.G       # activation energy for creep (J/mol)
        param["M"]  = r.M       # gas constant
    end
    return param
end

function make_param_saved_var()
    ps = ParamSaved()
    return Dict{Any,Any}(
        "maxTime"         => ps.maxTime,
        "lengthTime"      => ps.lengthTime,
        "switch_Tprofile" => ps.switch_Tprofile)
end

function make_sw()
        # eruption/cooling_module/viscous_relaxation control
        s = SW()
        return Dict{Any, Any}(
            "heat_cond"  => s.heat_cond,   # switch cooling module on/off
            "visc_relax" => s.visc_relax,   # switch viscous relaxation on/off
            "eruption"   => s.eruption
    )
end

function make_param_IC_Finder()
    # IC Finder parameters
    pic = ParamICFinder()
    param_IC_Finder = Dict{Any, Any}([])
    param_IC_Finder["max_count"]       = pic.max_count
    param_IC_Finder["min_eps_g"]       = pic.min_eps_g
    param_IC_Finder["eps_g_guess_ini"] = pic.eps_g_guess_ini
    param_IC_Finder["X_co2_guess_ini"] = pic.X_co2_guess_ini
    param_IC_Finder["fraction"]        = pic.fraction
    param_IC_Finder["delta_X_co2"]     = pic.delta_X_co2
    return param_IC_Finder
end