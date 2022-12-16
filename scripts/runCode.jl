using Chamber
include("./solver_methods.jl")

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
    log_vfr::Float64, depth::Number, methods::Dict=methods, method::String="Tsit5", odesetting=OdeSetting(), ini_eps_x::Float64=0.15, rheol::String="old") # ("silicic", 1e9, 0.2, 0.04, 0.001, -3.3, 8e3)
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

    to = get_timer("share")
    @timeit to "chamber" begin

    rc = rheol_composition_dict[composition]
    param =  Param{Float64}(composition=composition,
                            rheol=rheol,
                            rho_m0=rc.rho_m0,
                            rho_x0=rc.rho_x0,
                            c_m=rc.c_m,
                            c_x=rc.c_x,
                            L_m=rc.L_m)
    r = rheol_dict[rheol]
    if rheol == "new"
        param.A, param.B, param.G = r.A, r.B, r.G
    elseif rheol == "old"
        param.nn, param.AA, param.G, param.M = r.nn, r.AA, r.G, r.M
    end

    c = Const()
    param_IC_Finder = ParamICFinder()
    param_saved_var = ParamSaved()
    sw = SW()

    # Initial temperature and viscosity profile around chamber
    param_saved_var.storeSumk = zeros(param.maxn)
    param_saved_var.storeSumk_2 = zeros(param.maxn)
    param_saved_var.storeSumk_old = zeros(param.maxn)
    param_saved_var.storeSumk_2_old = zeros(param.maxn)

    volume_km3    = 10^log_volume_km3                 # range of volume in km3
    range_radius  = 1000*(volume_km3/(4*pi/3))^(1/3)  # range of radius in m
    V_0           = 4*pi/3*range_radius^3             # initial volume of the chamber (m^3)

    # thermal gradient
    param.Tb      = c.T_surface+c.T_gradient*depth        # background temperature crust (K)

    # lithostatic pressure
    P_0           = param.rho_r*c.grav_acc*depth   # initial chamber pressure (Pa)
    param.P_lit   = P_0
    param.P_lit_0 = P_0
    if param.single_eruption
        P_0 = P_0 + param.DP_crit
    end

    T_0 = find_liq(composition, InitialConc_H2O, InitialConc_CO2, P_0, ini_eps_x)

    T_in = T_0 + 50        # Temperature of inflowing magma (K)
    param.T_in = T_in

    # set the mass inflow rate
    param.Mdot_in_pass = build_mdot_in(param.fluxing, rc.rho_m0, log_vfr, P_0, T_in)

    rho_g0 = eos_g_rho_g(P_0, T_0)   # initial gas density
    eps_x0 = crystal_fraction_eps_x(composition, T_0, P_0, InitialConc_H2O, InitialConc_CO2)

    eps_m0 = 1-eps_x0
    rho = rc.rho_m0*eps_m0 + rc.rho_x0*eps_x0
    M_tot =  V_0*rho   # Total mass, initial
    M_co2_0 = InitialConc_CO2*V_0*rho   # Total mass of CO2, initial
    M_h2o_0 = InitialConc_H2O*V_0*rho   # Total mass of H2O, initial

    # IC Finder parameters
    param_IC_Finder.Tol = if composition == "silicic" 1e-9 else 1e-8 end

    if composition == "silicic"
        eps_g0, X_co20, C_co2, phase = IC_Finder_silicic(M_h2o_0, M_co2_0, M_tot, P_0, T_0, V_0, rc.rho_m0, param.mm_co2, param.mm_h2o, param_IC_Finder)
    elseif composition == "mafic"
        eps_g0, X_co20, C_co2, phase = IC_Finder_mafic(M_h2o_0, M_co2_0, M_tot, P_0, T_0, V_0, rc.rho_m0, param.mm_co2, param.mm_h2o, param_IC_Finder)
    end

    println("IC_Finder done: [eps_g0, X_co20, C_co2] = [$eps_g0, $X_co20, $C_co2]")
    println("phase: ", phase)
    @info("IC_Finder done: $eps_g0, $X_co20, $C_co2, $phase")
    param_saved_var.phase = phase

    # update initial bulk density (kg/m^3)
    rho_0  = rho_0_f(eps_g0, eps_x0, rho_g0, rc.rho_m0, rc.rho_x0)

    # update solubility
    if phase == 2
        X_co20 = 0.0
    end

    # Calculate the water content (concentration) for inflowing magma contents
    tot_h2o_frac_in = M_h2o_0/(rho_0*V_0)   # CHANGE MASS FRACTION FROM XCO2_IN
    tot_co2_frac_in = M_co2_0/(rho_0*V_0)
    if param.fluxing
        tot_co2_frac_in = param.XCO2_in*param.mm_co2/(param.XCO2_in*param.mm_co2+(1-param.XCO2_in)*param.mm_h2o)
        tot_h2o_frac_in = 1-tot_co2_frac_in
    end
    param.tot_h2o_frac_in = tot_h2o_frac_in
    param.tot_co2_frac_in = tot_co2_frac_in

    tot_Mass_0 = V_0*rho_0
    tot_Mass_H2O_0 = M_h2o_0
    tot_Mass_CO2_0 = M_co2_0

    if param.single_eruption
        sw.eruption = 1
    end

    # initialize vector to store quantities
    param_saved_var.storeTime = Vector{Float64}([0])
    param_saved_var.storeTemp = Vector{Float64}([T_0])

    # println(param)
    # println(param_IC_Finder)
    @info("sw: $(sw)")
    @info("IC_Finder parameters: $(param_IC_Finder)")
    @info("params: $(param)")


    """
        odeChamber(du,u,param,t)

    Define the ODE equation.
    """
    function odeChamber(du, u, param, t)
        composition = param.composition
        storeTime = param_saved_var.storeTime
        storeTemp = param_saved_var.storeTemp
        phase = param_saved_var.phase
        c_x, c_m = param.c_x, param.c_m
        L_e, L_m = param.L_e, param.L_m
        mm_h2o, mm_co2 = param.mm_h2o, param.mm_co2
        T_in = param.T_in
        P_lit_0, dP_lit_dt, dP_lit_dt_0, P_lit_drop_max = param.P_lit_0, param.dP_lit_dt, param.dP_lit_dt_0, param.P_lit_drop_max

        P0plusDP = u[1]
        T        = u[2]
        eps_g    = u[3]
        X_co2    = u[7]
        P_lit = P_lit_0 + dP_lit_dt_0*t
        if P_lit < P_lit_0-P_lit_drop_max
            P_lit = P_lit_0-P_lit_drop_max
            dP_lit_dt = 0.0
            param.dP_lit_dt = 0.0
        end
        P = P_lit + P0plusDP - P_lit_0
        # effective gas molar mass
        m_g = mm_co2*X_co2 + mm_h2o*(1-X_co2)

        if storeTime[end] == t
            storeTemp[end] = T
        elseif t != 0
            storeTime = [storeTime; t]
            storeTemp = [storeTemp; T]
            # push!(storeTime, t)
            # push!(storeTemp, T)
        end

        cross = findfirst(!=(0), diff(sign.(diff(storeTime))))
        if cross !== nothing
            cross_time= storeTime[end]
            storeTemp = [storeTemp[storeTime.<cross_time]; storeTemp[end]]
            storeTime = [storeTime[storeTime.<cross_time]; cross_time]
        end

        V, dV_dP, dV_dT             = compute_dXdP_dXdT(u[4], param, "r")
        rho_m, drho_m_dP, drho_m_dT = compute_dXdP_dXdT(u[5], param, "m")
        rho_x, drho_x_dP, drho_x_dT = compute_dXdP_dXdT(u[6], param, "x")
        eos_g_results = eos_g(P, T)
        rho_g, drho_g_dP, drho_g_dT = eos_g_results.rho_g, eos_g_results.drho_g_dP, eos_g_results.drho_g_dT

        total_Mass, M_h2o, M_co2    = u[8], u[9], u[10]
        m_h2o = M_h2o/(total_Mass)
        m_co2 = M_co2/(total_Mass)

        eps_x, deps_x_dP, deps_x_dT, deps_x_deps_g, deps_x_dmco2_t, deps_x_dmh2o_t =
            crystal_fraction(composition, T, P, m_h2o, m_co2)
        eps_m = 1-eps_x-eps_g

        m_eq, dm_eq_dP, dm_eq_dT, dm_eq_dX_co2, C_co2_t, dC_co2_dP, dC_co2_dT, dC_co2_dX_co2 =
            exsolve(composition, P, T, X_co2)

        if phase == 3
            C_co2 = C_co2_t
        elseif phase == 2
            C_co2 = m_co2
        end

        # specific heat of gas
        c_g = gas_heat_capacity(X_co2)

        rho, drho_dP, drho_dT, drho_deps_g, rc, drc_dP, drc_dT =
            build_rho_rc(eps_m, eps_g, eps_x, rho_m, rho_g, rho_x, drho_m_dP, drho_g_dP, drho_x_dP, 
            drho_m_dT, drho_g_dT, drho_x_dT, c_x, c_m, c_g, deps_x_dP, deps_x_dT)
        c = (1/rho)*(rho_x*eps_x*c_x+rho_m*eps_m*c_m+rho_g*eps_g*c_g)

        # boundary conditions
        Mdot_in, Mdot_out, Mdot_v_in, Mdot_v_out,Mdot_c_in, Mdot_c_out, Hdot_in, Hdot_out, P_loss, eta_r = 
            boundary_conditions_new(P, T, V, rho_m, rho_x, c, sw, T_in, M_h2o, M_co2, total_Mass, param, param_saved_var)

        A, b = build_matrix(phase, rho, drho_dP, V, dV_dP, drho_dT, dV_dT, drc_dP, rc, L_m, eps_x, drho_x_dP, T, 
            rho_x, deps_x_dP, L_e, dm_eq_dP, rho_m, eps_m, m_eq, drho_m_dP, drc_dT, drho_x_dT, deps_x_dT, dm_eq_dT, 
            drho_m_dT, Mdot_in, Mdot_out, P_loss, deps_x_dmh2o_t, m_h2o, m_co2, deps_x_dmco2_t, dP_lit_dt, 
            Hdot_in, Hdot_out, c_x, c_m, drho_deps_g, X_co2, m_g, eps_g, mm_h2o, drho_g_dP, rho_g, drho_g_dT, dm_eq_dX_co2, 
            mm_co2, c_g, dC_co2_dP, C_co2, dC_co2_dT, dC_co2_dX_co2, Mdot_v_in, Mdot_v_out, Mdot_c_in, Mdot_c_out)

        # coefficients in the system of unknowns Ax = B, here x= [dP/dt dT/dt deps_g/dt dX_co2/dt]
        # note: P, T, and phi are y(1), y(2) and y(3) respectively
        if phase == 3
            dDP_dt, dT_dt, deps_g_dt, dX_co2_dt = A\b
        elseif phase == 2
            dDP_dt, dT_dt = A\b
            deps_g_dt, dX_co2_dt = 0.0, 0.0
        end
        dP_dt          = dDP_dt + dP_lit_dt

        du[1]        = dDP_dt
        du[2]        = dT_dt
        du[3]        = deps_g_dt
        du[4]        = dV_dP*dP_dt + dV_dT*dT_dt + V*P_loss
        du[5]        = drho_m_dP*dP_dt + drho_m_dT*dT_dt
        du[6]        = drho_x_dP*dP_dt + drho_x_dT*dT_dt
        du[7]        = dX_co2_dt
        du[8]        = Mdot_in - Mdot_out
        du[9]        = Mdot_v_in - Mdot_v_out
        du[10]       = Mdot_c_in - Mdot_c_out
        param_saved_var.storeTime = storeTime
        param_saved_var.storeTemp = storeTemp
        return du
    end

    """
        stopChamber_MT(out,u,t,int)

    Define the stopping criteria for ODE solver. 
    """
    function stopChamber_MT(out, u, t, int)
        P_lit = param.P_lit
        DP_crit = param.DP_crit
        P0plusDP = u[1]
        T     = u[2]
        eps_g = u[3]
        V     = u[4]
        rho_m = u[5]
        tot_m = u[8]
        tot_w = u[9]
        tot_c = u[10]

        P = P0plusDP + P_lit - param.P_lit_0

        m_h20 = tot_w/tot_m
        m_co2 = tot_c/tot_m

        eps_x = crystal_fraction_eps_x(param.composition,T,P,m_h20,m_co2)
        m_eq_max = exsolve_meq(param.composition, P, T, 0.0)

        # MT's new stuff
        eps_m0 = 1 - eps_x

        m_h2o_melt = tot_w/(V*rho_m*eps_m0)
        m_co2_melt = tot_c/(V*rho_m*eps_m0)

        if param.composition == "silicic"
            C_co2_sat = exsolve3_silicic(P,T, m_h2o_melt)[1]
        elseif param.composition == "mafic"
            C_co2_sat = exsolve3_mafic(P,T, m_h2o_melt)[1]
        end

        out[1] = eps_x
        out[2] = eps_x/(1-eps_g)-0.8
        out[3] = if sw.eruption == 0 (P-P_lit)-DP_crit else -DP_crit end
        out[4] = if sw.eruption == 1 P_lit-P else -DP_crit end
        out[5] = eps_x-0.5
        out[6] = m_h2o_melt - m_eq_max
        out[7] = -(P0plusDP-param.P_lit_0+DP_crit)
        out[8] = m_co2_melt - C_co2_sat
    end

    """
        affect!(int, idx)

    Re-initialize the condition when the event happen.
    int, idx defined by DifferentialEquations 
    """
    function affect!(int, idx)
        println("*event idx: ", idx)
        storeTime = param_saved_var.storeTime
        storeTemp = param_saved_var.storeTemp
        storeTemp = storeTemp[storeTime.<int.t]
        storeTime = storeTime[storeTime.<int.t]
        param_saved_var.storeTime = storeTime
        param_saved_var.storeTemp = storeTemp

        if param.dP_lit_dt_0 == 0
            temp_P_lit = 0.0
        else
            if int.t <= abs(param.P_lit_drop_max/param.dP_lit_dt_0)
                temp_P_lit = param.dP_lit_dt_0*int.t
            else
                temp_P_lit = -param.P_lit_drop_max
            end
        end
        P_0 = int.u[1] + temp_P_lit

        m_h2o = int.u[9]/int.u[8]
        m_co2 = int.u[10]/int.u[8]

        eps_x0 =  crystal_fraction_eps_x(param.composition, int.u[2], P_0, m_h2o, m_co2)

        if idx == 3 && eps_x0 < 0.5
            sw.eruption = 1
            println("reached critical pressure and need to start an eruption,  time: ", int.t)
        elseif idx == 4
            sw.eruption = 0
            println("If it just finished an eruption...  time: ", int.t)
        elseif idx == 6 || idx == 8
            phase_here = param_saved_var.phase
            println("starting ic finder for conversion of phase,  time: $(int.t), phase_here: $phase_here")
            if param.composition == "silicic"
                eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder_silicic(int.u[9], int.u[10], int.u[8], P_0, int.u[2], int.u[4], int.u[5], param.mm_co2, param.mm_h2o, param_IC_Finder)
            elseif param.composition == "mafic"
                eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder_mafic(int.u[9], int.u[10], int.u[8], P_0, int.u[2], int.u[4], int.u[5], param.mm_co2, param.mm_h2o, param_IC_Finder)
            end

            param_saved_var.phase = phase
            if phase_here != phase
                println("1st try in IC Finder successful")
                int.u[3] = eps_g_temp
                int.u[7] = X_co2_temp
                C_co2 = C_co2_temp
            else
                println("trying new IC parameters...")
                param_IC_Finder.max_count = 150

                if param.composition == "silicic"
                    eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder_silicic(int.u[9], int.u[10], int.u[8], P_0, int.u[2], int.u[4], int.u[5], param.mm_co2, param.mm_h2o, param_IC_Finder)
                elseif param.composition == "mafic"
                    eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder_mafic(int.u[9], int.u[10], int.u[8], P_0, int.u[2], int.u[4], int.u[5], param.mm_co2, param.mm_h2o, param_IC_Finder)
                end
                param_saved_var.phase = phase
                ## change back to initial max_count
                param_IC_Finder.max_count = 100
                if phase_here != phase
                    println("2nd try in IC Finder successful")
                    int.u[3] = eps_g_temp
                    int.u[7] = X_co2_temp
                    C_co2 = C_co2_temp
                else
                    println("2nd try in IC Finder not successful, trying new IC parameters...")
                    param_IC_Finder.max_count = 100
                    param_IC_Finder.Tol = param_IC_Finder.Tol*0.1
                    if param.composition == "silicic"
                        eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder_silicic(int.u[9], int.u[10], int.u[8], P_0, int.u[2], int.u[4], int.u[5], param.mm_co2, param.mm_h2o, param_IC_Finder)
                    elseif param.composition == "mafic"
                        eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder_mafic(int.u[9], int.u[10], int.u[8], P_0, int.u[2], int.u[4], int.u[5], param.mm_co2, param.mm_h2o, param_IC_Finder)
                    end
                    param_saved_var.phase = phase
                    ## change back to initial Tol
                    param_IC_Finder.Tol = param_IC_Finder.Tol*10
                    if phase_here != phase
                        println("3rd try in IC Finder successful")
                        int.u[3] = eps_g_temp
                        int.u[7] = X_co2_temp
                        C_co2 = C_co2_temp
                    else
                        @warn("3rd try in IC Finder not successful")
                    end
                end
            end
            println("phase_here: ", phase_here, "  new_phase: ", phase)

        elseif idx == 1 || idx == 2 || idx == 5 || idx == 7 || idx === nothing
            if idx == 1
                println("eps_x became 0.")
            elseif idx == 2
                println("eps_x/(1-eps_g) became 0.8")
            elseif idx == 5
                println("eps_x became 0.5")
            elseif idx == 7
                println("too much underpressure - collapse")
            elseif idx === nothing
                println("you reached the end of time")
            end
            terminate!(int)
        end
    end

    tspan    = (0, end_time)
    IC       = [P_0, T_0, eps_g0, V_0, rc.rho_m0, rc.rho_x0, X_co20, tot_Mass_0, tot_Mass_H2O_0, tot_Mass_CO2_0]
    println("tspan: ", tspan)
    println("IC: ", IC)
    @info("IC: $IC")
    cb = VectorContinuousCallback(stopChamber_MT, affect!, 8, rootfind=SciMLBase.RightRootFind)
    prob = ODEProblem(odeChamber,IC,tspan,param)
    sol = solve(prob, methods[method], callback=cb, reltol=odesetting.reltol, abstol=odesetting.abstol, dt=odesetting.first_step, dtmax=odesetting.max_step)
    end
    println(to)
    @info(to)
    close(io)
    write_csv(sol, path)
    plot_figs("$path/out.csv", path)

    println(".. Done!")
    return sol
end

# chamber("silicic", 1e9, 0.2, 0.04, 0.001, -3.3, 8e3)
