# Use for runCode.jl, solving the ode problems. 
"""
    odeChamber(du,u,param,t)

Define the ODE equation.
"""
function odeChamber(du, u, param, t)
    storeTime = param_saved_var["storeTime"]
    storeTemp = param_saved_var["storeTemp"]
    phase = param_saved_var["phase"]

    P0plusDP = u[1]
    T        = u[2]
    eps_g    = u[3]
    X_co2    = u[7]
    P_lit = param["P_lit_0"]+param["dP_lit_dt_0"]*t
    if P_lit < param["P_lit_0"]-param["P_lit_drop_max"]
        P_lit = param["P_lit_0"]-param["P_lit_drop_max"]
        param["dP_lit_dt"] = 0
    end
    P = P_lit + P0plusDP - param["P_lit_0"]
    # effective gas molar mass
    m_g = param["mm_co2"]*X_co2+param["mm_h2o"]*(1-X_co2)

    if ~isempty(storeTime)
        if storeTime[end] == t
            storeTime[end] = t
            storeTemp[end] = T
        elseif t == 0
        else
            storeTime = [storeTime; t]
            storeTemp = [storeTemp; T]
        end
    elseif t == 0
    else
        storeTime = [storeTime; t]
        storeTemp = [storeTemp; T]
    end
    cross = findfirst(>(0),abs.(diff(sign.(diff(storeTime)))))
    if cross !== nothing
        cross_time=  storeTime[end]
        storeTemp = [storeTemp[storeTime.<cross_time]; storeTemp[end]]
        storeTime = [storeTime[storeTime.<cross_time]; cross_time]
    end

    V              = u[4]
    dV_dP          = V/param["beta_r"]
    dV_dT          = -V*param["alpha_r"]

    rho_m          = u[5]
    drho_m_dP      = rho_m/param["beta_m"]
    drho_m_dT      = -rho_m*param["alpha_m"]

    rho_x          = u[6]
    drho_x_dP      = rho_x/param["beta_x"]
    drho_x_dT      = -rho_x*param["alpha_x"]

    eos_g_results = eos_g(P,T)
    rho_g,drho_g_dP,drho_g_dT = eos_g_results["rho_g"],eos_g_results["drho_g_dP"],eos_g_results["drho_g_dT"]

    M_h2o    = u[9]
    M_co2    = u[10]
    total_Mass = u[8]

    m_h2o= M_h2o/(total_Mass)
    m_co2=M_co2/(total_Mass)

    eps_x, deps_x_dP, deps_x_dT, deps_x_deps_g, deps_x_dmco2_t, deps_x_dmh2o_t = crystal_fraction(param["composition"],T,P,m_h2o,m_co2)

    eps_m=1-eps_x-eps_g
    rho = eps_m*rho_m + eps_g*rho_g + eps_x*rho_x

    if param["composition"] == "silicic"
        m_eq,dm_eq_dP,dm_eq_dT,dm_eq_dX_co2, C_co2_t,dC_co2_dP, dC_co2_dT, dC_co2_dX_co2 = exsolve_silicic(P,T, X_co2)
    elseif param["composition"] == "mafic"
        m_eq,dm_eq_dP,dm_eq_dT,dm_eq_dX_co2, C_co2_t,dC_co2_dP, dC_co2_dT, dC_co2_dX_co2 = exsolve_mafic(P,T, X_co2)
    end

    if phase == 3
        C_co2=C_co2_t
    elseif phase == 2
        C_co2 = m_co2
    end

    rho         = eps_m*rho_m + eps_g*rho_g + eps_x*rho_x
    drho_dP     = eps_m*drho_m_dP + eps_g*drho_g_dP + eps_x*drho_x_dP+(rho_x-rho_m)*deps_x_dP
    drho_dT     = eps_m*drho_m_dT + eps_g*drho_g_dT + eps_x*drho_x_dT+(rho_x-rho_m)*deps_x_dT
    drho_deps_g = -rho_m + rho_g

    # % specific heat of gas
    c_g = gas_heat_capacity(X_co2)[1]
    c = (1/rho)*(rho_x*eps_x*param["c_x"]+rho_m*eps_m*param["c_m"]+rho_g*eps_g*c_g)

    # computing the product of density and specific heat for the mixture and
    # its derivatives
    rc              = rho_x*eps_x*param["c_x"]+rho_m*eps_m*param["c_m"]+rho_g*eps_g*c_g
    drc_dP          = eps_x*param["c_x"]*drho_x_dP+eps_g*c_g*drho_g_dP+eps_m*param["c_m"]*drho_m_dP+(rho_x*param["c_x"]-rho_m*param["c_m"])*deps_x_dP
    drc_dT          = eps_x*param["c_x"]*drho_x_dT+eps_g*c_g*drho_g_dT+eps_m*param["c_m"]*drho_m_dT+(rho_x*param["c_x"]-rho_m*param["c_m"])*deps_x_dT

    # boundary conditions
    Mdot_in, Mdot_out, Mdot_v_in, Mdot_v_out,Mdot_c_in, Mdot_c_out, Hdot_in, Hdot_out, P_loss,eta_r = boundary_conditions_new(P, T, V, rho_m, rho_x, c, sw, param["T_in"], M_h2o, M_co2, total_Mass, param, param_saved_var)
    # coefficients in the system of unknowns Ax = B, here x= [dP/dt dT/dt deps_g/dt dX_co2/dt]
    # note: P, T, and phi are y(1), y(2) and y(3) respectively
    # values matrix A
    # conservation of (total) mass
    a11 = (1/rho)*drho_dP +(1/V)*dV_dP
    a12 = (1/rho)*drho_dT + (1/V)*dV_dT
    a13 = (1/rho)*drho_deps_g
    a14 = 0

    # conservation of water mass
    a21 = eps_m*dm_eq_dP-m_eq*deps_x_dP+m_eq*eps_m*dV_dP/V+m_eq*eps_m* drho_m_dP/rho_m+(1-X_co2)/m_g*eps_g*param["mm_h2o"]*drho_g_dP/rho_m+(1-X_co2)/m_g*eps_g*rho_g*param["mm_h2o"]*dV_dP/(rho_m*V)
    a22 = eps_m*dm_eq_dT-m_eq*deps_x_dT-m_eq*eps_m*(1/V)*dV_dT+m_eq*eps_m*drho_m_dT/rho_m+(1-X_co2)/m_g*eps_g*param["mm_h2o"]*drho_g_dT/rho_m+(1-X_co2)/m_g*eps_g*rho_g*param["mm_h2o"]*dV_dT/(rho_m*V)
    a23 = -m_eq+(1-X_co2)/m_g*rho_g*param["mm_h2o"]/rho_m
    a24 = real(eps_m*dm_eq_dX_co2-1/m_g*eps_g*rho_g*param["mm_h2o"]/rho_m-(1-X_co2)*eps_g*rho_g*param["mm_h2o"]*(param["mm_co2"]-param["mm_h2o"])/(Complex(m_g)^2*rho_m))

    # conservation of (total) enthalpy - all divided by rc*T*V
    a31  = drc_dP/(rc)+dV_dP/V+param["L_m"]*(-eps_x*drho_x_dP/(rc*T)-rho_x*deps_x_dP/(rc*T)-rho_x*eps_x*dV_dP/(rc*V*T))+param["L_e"]*(-dm_eq_dP*rho_m*eps_m/(rc*T)-m_eq*eps_m*drho_m_dP/(rc*T)+m_eq*rho_m*deps_x_dP/(rc*T)-m_eq*rho_m*eps_m*dV_dP/(rc*V*T))
    a32  = drc_dT/(rc)+1/T+dV_dT/V+param["L_m"]*(-eps_x*drho_x_dT/(rc*T)-rho_x*deps_x_dT/(rc*T)-rho_x*eps_x*dV_dT/(rc*T*V))+param["L_e"]*(-dm_eq_dT*rho_m*eps_m/(rc*T)-m_eq*eps_m*drho_m_dT/(rc*T)+m_eq*rho_m*deps_x_dT/(rc*T)-m_eq*rho_m*eps_m*dV_dT/(rc*T*V))
    a33  = (rho_g*c_g-rho_m*param["c_m"])/rc +param["L_e"]*m_eq*rho_m/(rc*T)
    a34  = -param["L_e"]*rho_m*eps_m*dm_eq_dX_co2/(rc*T)

    # conservation of CO2 mass
    a41 = eps_m*dC_co2_dP-C_co2*deps_x_dP+C_co2*eps_m*dV_dP/V+C_co2*eps_m*drho_m_dP/rho_m+X_co2/m_g*eps_g*param["mm_co2"]*drho_g_dP/rho_m+X_co2/m_g*eps_g*rho_g*param["mm_co2"]*dV_dP/(rho_m*V)
    a42 = eps_m*dC_co2_dT-C_co2*deps_x_dT-C_co2*eps_m*(1/V)*dV_dT+C_co2*eps_m*drho_m_dT/rho_m+X_co2/m_g*eps_g*param["mm_co2"]*drho_g_dT/rho_m+X_co2/m_g*eps_g*rho_g*param["mm_co2"]*dV_dT/(rho_m*V)
    a43 = -C_co2+X_co2/m_g*rho_g*param["mm_co2"]/rho_m
    a44 = real(eps_m*dC_co2_dX_co2+1/m_g*eps_g*rho_g*param["mm_co2"]/rho_m-X_co2*eps_g*rho_g*param["mm_co2"]*(param["mm_co2"]-param["mm_h2o"])/(Complex(m_g)^2*rho_m))

    # values vector B
    # conservation of (total) mass
    dM_h2o_t_dt=1/(rho*V)*((Mdot_v_in-Mdot_v_out)-m_h2o*(Mdot_in-Mdot_out))
    dM_co2_t_dt=1/(rho*V)*((Mdot_c_in-Mdot_c_out)-m_co2*(Mdot_in-Mdot_out))

    b1  =  (Mdot_in - Mdot_out)/(rho*V) - P_loss-(rho_x-rho_m)/rho*deps_x_dmh2o_t*dM_h2o_t_dt-(rho_x-rho_m)/rho*deps_x_dmco2_t*dM_co2_t_dt-a11*param["dP_lit_dt"]
    # conservation of water mass
    b2  = (Mdot_v_in - Mdot_v_out)/(rho_m*V)-m_eq*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)-m_eq*eps_m*P_loss-(1-X_co2)/m_g*rho_g/rho_m*eps_g*param["mm_h2o"]*P_loss-a21*param["dP_lit_dt"]
    # conservation of (total) enthalpy
    b3  =  (Hdot_in - Hdot_out)/(rc*T*V) - 1/(rc*V*T)*((rho_x*param["c_x"]-rho_m*param["c_m"])*T*V*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt))+param["L_m"]*rho_x*V/(rc*V*T)*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)+param["L_e"]*m_eq*rho_m*V/(rc*V*T)*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)+P_loss*(-1+param["L_m"]*rho_x*eps_x*V/(rc*V*T)+param["L_e"]*m_eq*rho_m*eps_m*V/(rc*V*T))-a31*param["dP_lit_dt"]
    # conservation of CO2 mass
    b4  =  (Mdot_c_in - Mdot_c_out)/(rho_m*V)-C_co2*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)-C_co2*eps_m*P_loss-(X_co2)/m_g*rho_g/rho_m*eps_g*param["mm_co2"]*P_loss-a41*param["dP_lit_dt"]

    if phase == 3
        # set up matrices to solve using Cramer"s rule
        A          = [a11 a12 a13 a14; a21 a22 a23 a24; a31 a32 a33 a34; a41 a42 a43 a44]
        A_P        = [b1  a12 a13 a14; b2  a22 a23 a24; b3  a32 a33 a34; b4  a42 a43 a44]
        A_T        = [a11 b1  a13 a14; a21 b2  a23 a24; a31 b3  a33 a34; a41 b4  a43 a44]
        A_eps_g    = [a11 a12 b1  a14; a21 a22 b2  a24; a31 a32 b3  a34; a41 a42 b4  a44]
        A_X_co2    = [a11 a12 a13 b1 ; a21 a22 a23 b2 ; a31 a32 a33 b3 ; a41 a42 a43 b4 ]

        det_A          = det(A)
        dDP_dt         = det(A_P)/det_A
        dT_dt          = det(A_T)/det_A
        deps_g_dt      = det(A_eps_g)/det_A
        dX_co2_dt      = det(A_X_co2)/det_A
        dP_dt          = dDP_dt+param["dP_lit_dt"]
        dV_dt          = dV_dP*dP_dt + dV_dT*dT_dt + V*P_loss
        drho_m_dt      = drho_m_dP*dP_dt + drho_m_dT*dT_dt
        drho_x_dt      = drho_x_dP*dP_dt + drho_x_dT*dT_dt
    else
        A          = [a11 a12; a31 a32]
        A_P        = [b1  a12; b3  a32]
        A_T        = [a11 b1 ; a31 b3]
        
        det_A          = det(A)
        dDP_dt         = det(A_P)/det_A
        dT_dt          = det(A_T)/det_A
        deps_g_dt      = 0
            
        dP_dt          = dDP_dt+param["dP_lit_dt"]
        dV_dt          = dV_dP*dP_dt + dV_dT*dT_dt + V*P_loss
        drho_m_dt      = drho_m_dP*dP_dt + drho_m_dT*dT_dt
        drho_x_dt      = drho_x_dP*dP_dt + drho_x_dT*dT_dt
        dX_co2_dt      = 0
    end
    # write(io, "  $t,  $(u[1]),  $(u[3]),  $(u[7])\n")
    du[1]        = dDP_dt
    du[2]        = dT_dt
    du[3]        = deps_g_dt
    du[4]        = dV_dt
    du[5]        = drho_m_dt
    du[6]        = drho_x_dt
    du[7]        = dX_co2_dt
    du[8]        = Mdot_in - Mdot_out
    du[9]        = Mdot_v_in - Mdot_v_out
    du[10]       = Mdot_c_in - Mdot_c_out
    param_saved_var["storeTime"] = storeTime
    param_saved_var["storeTemp"] = storeTemp
    return du
end

"""
    stopChamber_MT(out,u,t,int)

Define the stopping criteria for ODE solver. 
"""
function stopChamber_MT(out, u, t, int)
    P_lit = param["P_lit"]
    DP_crit = param["DP_crit"]
    P0plusDP = u[1]
    T     = u[2]
    eps_g = u[3]
    V     = u[4]
    rho_m = u[5]
    tot_m = u[8]
    tot_w = u[9]
    tot_c = u[10]

    P = P0plusDP + P_lit - param["P_lit_0"] 

    m_h20 = tot_w/tot_m
    m_co2 = tot_c/tot_m

    eps_x = crystal_fraction_eps_x(param["composition"],T,P,m_h20,m_co2)
    if param["composition"] == "silicic"
        m_eq_max = exsolve_silicic(P, T, 0)[1]

    elseif param["composition"] == "mafic"
        m_eq_max = exsolve_mafic(P, T, 0)[1]
    end
    # MT's new stuff
    eps_m0 = 1 - eps_x

    m_h2o_melt = tot_w/(V*rho_m*eps_m0)
    m_co2_melt = tot_c/(V*rho_m*eps_m0)

    if param["composition"] == "silicic"
        C_co2_sat = exsolve3_silicic(P,T, m_h2o_melt)[1]
    elseif param["composition"] == "mafic"
        C_co2_sat = exsolve3_mafic(P,T, m_h2o_melt)[1]
    end

    out[1] = eps_x
    out[2] = eps_x/(1-eps_g)-0.8
    out[3] = if sw.eruption == 0 (P-P_lit)-DP_crit else -DP_crit end
    out[4] = if sw.eruption == 1 P_lit-P else -DP_crit end
    out[5] = eps_x-0.5
    out[6] = m_h2o_melt - m_eq_max
    out[7] = -(P0plusDP-param["P_lit_0"]+DP_crit)
    out[8] = m_co2_melt - C_co2_sat
    param["out"] = out
end

"""
    affect!(int, idx)

Re-initialize the condition when the event happen.
int, idx defined by DifferentialEquations 
"""
function affect!(int, idx)
    println("*event idx: ", idx)
    # write(io, "*event idx: $idx \n")

    storeTime = param_saved_var["storeTime"]
    storeTemp = param_saved_var["storeTemp"]
    storeTemp = storeTemp[storeTime.<int.t]
    storeTime = storeTime[storeTime.<int.t]
    param_saved_var["storeTime"] = storeTime
    param_saved_var["storeTemp"] = storeTemp

    if param["dP_lit_dt_0"] == 0
        temp_P_lit = 0
    else
        if int.t <= abs(param["P_lit_drop_max"]/param["dP_lit_dt_0"])
            temp_P_lit = param["dP_lit_dt_0"]*int.t
        else
            temp_P_lit = -param["P_lit_drop_max"]
        end
    end
    P_0 = int.u[1] + temp_P_lit

    m_h2o = int.u[9]/int.u[8]
    m_co2 = int.u[10]/int.u[8]

    eps_x0 =  crystal_fraction_eps_x(param["composition"], int.u[2], P_0, m_h2o, m_co2)

    if idx == 3 && eps_x0 < 0.5
        sw.eruption = 1
        println("reached critical pressure and need to start an eruption,  time: ", int.t)
        if "out" in keys(param)
            # write(io, " stopChamber_MT: $(param["out"])\n")
        end
    elseif idx == 4
        sw.eruption = 0
        println("If it just finished an eruption...  time: ", int.t)
        if "out" in keys(param)
            # write(io, " stopChamber_MT: $(param["out"])\n")
        end
    elseif idx == 6 || idx == 8
        phase_here = param_saved_var["phase"]
        println("starting ic finder for conversion of phase,  time: $(int.t), phase_here: $phase_here")
        # write(io, "starting ic finder for conversion of phase,  time: $(int.t), phase_here: $phase_here\n")
        if "out" in keys(param)
            # write(io, " stopChamber_MT: $(param["out"])\n")
        end
        if param["composition"] == "silicic"
            eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder_silicic(int.u[9], int.u[10], int.u[8], P_0, int.u[2], int.u[4], int.u[5], param["mm_co2"], param["mm_h2o"], param_IC_Finder)
            # write(io, " 1. [$eps_g_temp, $X_co2_temp, $C_co2_temp, $phase] = IC_Finder_silicic($(int.u[9]), $(int.u[10]), $(int.u[8]), $P_0, $(int.u[2]), $(int.u[4]), $(int.u[5])),  max_count: $(param_IC_Finder["max_count"])\n")
        elseif param["composition"] == "mafic"
            eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder_mafic(int.u[9], int.u[10], int.u[8], P_0, int.u[2], int.u[4], int.u[5], param["mm_co2"], param["mm_h2o"], param_IC_Finder)
            # write(io, " 1. [$eps_g_temp, $X_co2_temp, $C_co2_temp, $phase] = IC_Finder_mafic($(int.u[9]), $(int.u[10]), $(int.u[8]), $P_0, $(int.u[2]), $(int.u[4]), $(int.u[5])),  max_count: $(param_IC_Finder["max_count"])\n")
        end

        param_saved_var["phase"] = phase
        if phase_here != phase
            println("1st try in IC Finder successful")
            # write(io, "1st try in IC Finder successful\n")
            int.u[3] = eps_g_temp
            int.u[7] = X_co2_temp
            C_co2 = C_co2_temp
        else
            println("trying new IC parameters...")
            # write(io, "trying new IC parameters...\n")
            param_IC_Finder["max_count"] = 150
            param_IC_Finder["Tol"] = if composition == "Silicic" 1e-9 else 1e-8 end
            param_IC_Finder["min_eps_g"] = 1e-10
            param_IC_Finder["eps_g_guess_ini"] = 1e-2
            param_IC_Finder["X_co2_guess_ini"] = 0.2
            param_IC_Finder["fraction"] = 0.2
            param_IC_Finder["delta_X_co2"] = 1e-2
            if param["composition"] == "silicic"
                eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder_silicic(int.u[9], int.u[10], int.u[8], P_0, int.u[2], int.u[4], int.u[5], param["mm_co2"], param["mm_h2o"], param_IC_Finder)
                # write(io, " 2. [$eps_g_temp, $X_co2_temp, $C_co2_temp, $phase] = IC_Finder_silicic($(int.u[9]), $(int.u[10]), $(int.u[8]), $P_0, $(int.u[2]), $(int.u[4]), $(int.u[5])),  max_count: $(param_IC_Finder["max_count"])\n")
            elseif param["composition"] == "mafic"
                eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder_mafic(int.u[9], int.u[10], int.u[8], P_0, int.u[2], int.u[4], int.u[5], param["mm_co2"], param["mm_h2o"], param_IC_Finder)
                # write(io, " 2. [$eps_g_temp, $X_co2_temp, $C_co2_temp, $phase] = IC_Finder_mafic($(int.u[9]), $(int.u[10]), $(int.u[8]), $P_0, $(int.u[2]), $(int.u[4]), $(int.u[5])),  max_count: $(param_IC_Finder["max_count"])\n")
            end
            param_saved_var["phase"] = phase
            ## change back to initial max_count
            param_IC_Finder["max_count"] = 100
            if phase_here != phase
                println("2nd try in IC Finder successful")
                # write(io, "2nd try in IC Finder successful\n")
                int.u[3] = eps_g_temp
                int.u[7] = X_co2_temp
                C_co2 = C_co2_temp
            else
                println("2nd try in IC Finder not successful, trying new IC parameters...")
                # write(io, "2nd try in IC Finder not successful, trying new IC parameters...\n")
                param_IC_Finder["max_count"] = 100
                param_IC_Finder["Tol"] = param_IC_Finder["Tol"]*0.1
                if param["composition"] == "silicic"
                    eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder_silicic(int.u[9], int.u[10], int.u[8], P_0, int.u[2], int.u[4], int.u[5], param["mm_co2"], param["mm_h2o"], param_IC_Finder)
                    # write(io, " 3. [$eps_g_temp, $X_co2_temp, $C_co2_temp, $phase] = IC_Finder_silicic($(int.u[9]), $(int.u[10]), $(int.u[8]), $P_0, $(int.u[2]), $(int.u[4]), $(int.u[5])),  max_count: $(param_IC_Finder["max_count"])\n")
                elseif param["composition"] == "mafic"
                    eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder_mafic(int.u[9], int.u[10], int.u[8], P_0, int.u[2], int.u[4], int.u[5], param["mm_co2"], param["mm_h2o"], param_IC_Finder)
                    # write(io, " 3. [$eps_g_temp, $X_co2_temp, $C_co2_temp, $phase] = IC_Finder_mafic($(int.u[9]), $(int.u[10]), $(int.u[8]), $P_0, $(int.u[2]), $(int.u[4]), $(int.u[5])),  max_count: $(param_IC_Finder["max_count"])\n")
                end
                param_saved_var["phase"] = phase
                ## change back to initial Tol
                param_IC_Finder["Tol"] = param_IC_Finder["Tol"]*10
                if phase_here != phase
                    println("3rd try in IC Finder successful")
                    # write(io, "3rd try in IC Finder successful\n")
                    int.u[3] = eps_g_temp
                    int.u[7] = X_co2_temp
                    C_co2 = C_co2_temp
                else
                    println("3rd try in IC Finder not successful")
                    # write(io, "3rd try in IC Finder not successful\n")
                end
            end
        end
        println("phase_here: ", phase_here, "  new_phase: ", phase)
        # write(io, " phase_here: $phase_here, new_phase: $phase\n")

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
    # write(io, "time: $(int.t)\n")
    # write(io, "IC = [$(int.u[1]), $(int.u[2]), $(int.u[3]), $(int.u[4]), $(int.u[5]), $(int.u[6]), $(int.u[7]), $(int.u[8]), $(int.u[9]), $(int.u[10])]\n")
end