# For runCode.jl, solving the ode problems. 
"""
    odeChamber(du::Vector{Float64}, u::Vector{Float64}, params::Tuple{Param{Float64}, ParamSaved{Float64}, SW{Int8}}, t::Float64)

- Define the ODE equation.
- Solve the model for eruption frequency of upper crustal silicic magma chambers using an ODE solver.

# Arguments
- `du`: An array that stores the output of the ODE solver, i.e., the values of the derivatives of the solution `u` with respect to time `t`.
- `u`: An array that stores the values of the solution at each time step `t`.
- `p`: A tuple that stores the model parameters and some saved variables, which are described in more detail below.
- `t`: The time points corresponding to the saved values of the ODE solution.

The arguments `du`, `u`, `p`, and `t` are from the DifferentialEquations.jl package. These argument formats are specific to the DifferentialEquations.jl package.

# Parameters
- `param`: A custom parameter containing physical constants and other model parameters.
- `param_saved_var`: A custom parameter used to store values from the previous time step.
- `sw`: A custom parameter used to control simulation behavior.

# Returns
The function modifies `du` in place to store the values of the derivatives of the solution `u` with respect to time `t`.
"""
function odeChamber(
    du::Vector{Float64},
    u::Vector{Float64},
    p::Tuple{Param{Float64},ParamSaved{Float64},SW{Int8}},
    t::Float64,
)
    param, param_saved_var, sw = p
    composition = param.composition
    storeTime = param_saved_var.storeTime
    storeTemp = param_saved_var.storeTemp
    phase = param_saved_var.phase
    c_x, c_m = param.c_x, param.c_m
    L_e, L_m = param.L_e, param.L_m
    mm_h2o, mm_co2 = param.mm_h2o, param.mm_co2
    T_in = param.T_in
    P_lit_0, dP_lit_dt, dP_lit_dt_0, P_lit_drop_max = param.P_lit_0,
    param.dP_lit_dt, param.dP_lit_dt_0,
    param.P_lit_drop_max

    P0plusDP = u[1]
    T = u[2]
    eps_g = u[3]
    X_co2 = u[7]
    P_lit = P_lit_0 + dP_lit_dt_0 * t
    if P_lit < P_lit_0 - P_lit_drop_max
        P_lit = P_lit_0 - P_lit_drop_max
        dP_lit_dt = 0.0
        param.dP_lit_dt = 0.0
    end
    P = P_lit + P0plusDP - P_lit_0
    # effective gas molar mass
    m_g = mm_co2 * X_co2 + mm_h2o * (1 - X_co2)

    #=
    NOTE:
        push! method is NOT work for odechamber, it may cause some errors when running `chamber`. 
        using `storeTime = [storeTime; t]` instead of `push!(storeTime, t)`
    =#
    if storeTime[end] == t
        storeTemp[end] = T
    elseif t != 0
        storeTime = [storeTime; t]
        storeTemp = [storeTemp; T]
    end

    cross = findfirst(!=(0), diff(sign.(diff(storeTime))))
    if cross !== nothing
        cross_time = storeTime[end]
        storeTemp = [storeTemp[storeTime .< cross_time]; storeTemp[end]]
        storeTime = [storeTime[storeTime .< cross_time]; cross_time]
    end

    V, dV_dP, dV_dT = compute_dXdP_dXdT(u[4], param, "r")
    rho_m, drho_m_dP, drho_m_dT = compute_dXdP_dXdT(u[5], param, "m")
    rho_x, drho_x_dP, drho_x_dT = compute_dXdP_dXdT(u[6], param, "x")
    eos_g_results = eos_g(P, T)
    rho_g, drho_g_dP, drho_g_dT = eos_g_results.rho_g,
    eos_g_results.drho_g_dP,
    eos_g_results.drho_g_dT

    total_Mass, M_h2o, M_co2 = u[8], u[9], u[10]
    m_h2o = M_h2o / (total_Mass)
    m_co2 = M_co2 / (total_Mass)

    eps_x, deps_x_dP, deps_x_dT, deps_x_deps_g, deps_x_dmco2_t, deps_x_dmh2o_t = crystal_fraction(
        composition, T, P, m_h2o, m_co2
    )
    eps_m = 1 - eps_x - eps_g

    m_eq, dm_eq_dP, dm_eq_dT, dm_eq_dX_co2, C_co2_t, dC_co2_dP, dC_co2_dT, dC_co2_dX_co2 = exsolve(
        composition, P, T, X_co2
    )

    if phase == 3
        C_co2 = C_co2_t
    elseif phase == 2
        C_co2 = m_co2
    end

    # specific heat of gas
    c_g = gas_heat_capacity(X_co2)

    rho, drho_dP, drho_dT, drho_deps_g, rc, drc_dP, drc_dT = build_rho_rc(
        eps_m,
        eps_g,
        eps_x,
        rho_m,
        rho_g,
        rho_x,
        drho_m_dP,
        drho_g_dP,
        drho_x_dP,
        drho_m_dT,
        drho_g_dT,
        drho_x_dT,
        c_x,
        c_m,
        c_g,
        deps_x_dP,
        deps_x_dT,
    )
    c = (1 / rho) * (rho_x * eps_x * c_x + rho_m * eps_m * c_m + rho_g * eps_g * c_g)

    # boundary conditions
    Mdot_in, Mdot_out, Mdot_v_in, Mdot_v_out, Mdot_c_in, Mdot_c_out, Hdot_in, Hdot_out, P_loss, eta_r = boundary_conditions_new(
        P, T, V, rho_m, rho_x, c, sw, T_in, M_h2o, M_co2, total_Mass, param, param_saved_var
    )

    A, b = build_matrix(
        phase,
        rho,
        drho_dP,
        V,
        dV_dP,
        drho_dT,
        dV_dT,
        drc_dP,
        rc,
        L_m,
        eps_x,
        drho_x_dP,
        T,
        rho_x,
        deps_x_dP,
        L_e,
        dm_eq_dP,
        rho_m,
        eps_m,
        m_eq,
        drho_m_dP,
        drc_dT,
        drho_x_dT,
        deps_x_dT,
        dm_eq_dT,
        drho_m_dT,
        Mdot_in,
        Mdot_out,
        P_loss,
        deps_x_dmh2o_t,
        m_h2o,
        m_co2,
        deps_x_dmco2_t,
        dP_lit_dt,
        Hdot_in,
        Hdot_out,
        c_x,
        c_m,
        drho_deps_g,
        X_co2,
        m_g,
        eps_g,
        mm_h2o,
        drho_g_dP,
        rho_g,
        drho_g_dT,
        dm_eq_dX_co2,
        mm_co2,
        c_g,
        dC_co2_dP,
        C_co2,
        dC_co2_dT,
        dC_co2_dX_co2,
        Mdot_v_in,
        Mdot_v_out,
        Mdot_c_in,
        Mdot_c_out,
    )

    # coefficients in the system of unknowns Ax = B, here x= [dP/dt dT/dt deps_g/dt dX_co2/dt]
    if phase == 3
        dDP_dt, dT_dt, deps_g_dt, dX_co2_dt = A \ b
    elseif phase == 2
        dDP_dt, dT_dt = A \ b
        deps_g_dt, dX_co2_dt = 0.0, 0.0
    end
    dP_dt = dDP_dt + dP_lit_dt

    du[1] = dDP_dt
    du[2] = dT_dt
    du[3] = deps_g_dt
    du[4] = dV_dP * dP_dt + dV_dT * dT_dt + V * P_loss
    du[5] = drho_m_dP * dP_dt + drho_m_dT * dT_dt
    du[6] = drho_x_dP * dP_dt + drho_x_dT * dT_dt
    du[7] = dX_co2_dt
    du[8] = Mdot_in - Mdot_out
    du[9] = Mdot_v_in - Mdot_v_out
    du[10] = Mdot_c_in - Mdot_c_out
    param_saved_var.storeTime = storeTime
    param_saved_var.storeTemp = storeTemp
    return du
end

"""
    stopChamber_MT(out, u::Vector{Float64}, t::Float64, int, sw::SW{Int8}, param::Param{Float64})

Define the stopping criteria for an ODE solver that simulates a magma chamber.

# Arguments:
- `out`: An array where the function should save the condition value at the right index. The maximum index of `out` should be specified in the `len` property of `callback`, which allows for a chain of `len` events, triggering the `ith` event when `out[i] = 0`. The function returns the value of `out[8]` as the last condition. Checking Event Handling and Callback Functions page of `DifferentialEquations.jl` for more details.
- `u`: A vector containing the state of the system at time `t`.
- `t`: The current time of the ODE solver.
- `int`: The current state of the integrator. It's format is from the DifferentialEquations.jl package
- `sw`: A custom parameter used to control simulation behavior.
- `param`: A custom parameter containing physical constants and other model parameters.

The arguments `out`, `u`, `t`, and `int` are from the DifferentialEquations.jl package. These argument formats are specific to the DifferentialEquations.jl package.

# Returns
The `out` array is modified in-place to contain the condition values at the current state of the system. The function computes several quantities based on the current state of the system and the model parameters, such as the crystal fraction, the maximum amount of exsolved fluid, and the amount of CO2 in the melt. It then uses these values to calculate the condition values to be stored in the `out` array, as follows:
- `out[1]`: the crystal fraction.
- `out[2]`: the ratio of crystal fraction to liquid fraction, minus 0.8.
- `out[3]`: if an eruption is not occurring, the pressure difference between the lithostatic pressure and the current pressure, minus the critical pressure difference. Otherwise, negative the critical pressure difference.
- `out[4]`: if an eruption is occurring, the lithostatic pressure minus the current pressure. Otherwise, negative the critical pressure difference.
- `out[5]`: the crystal fraction minus 0.5.
- `out[6]`: the difference between the amount of water in the melt and the maximum amount of exsolved fluid.
- `out[7]`: negative the difference between the initial lithostatic pressure and the current pressure plus the critical pressure difference.
- `out[8]`: the difference between the amount of CO2 in the melt and the saturation concentration.

Note that the `out` and `u` arguments are in the format expected by the DifferentialEquations.jl package, and the function is intended to be used as a condition for a callback function.
"""
function stopChamber_MT(
    out, u::Vector{Float64}, t::Float64, int, sw::SW{Int8}, param::Param{Float64}
)
    composition = param.composition
    P_lit = param.P_lit
    DP_crit = param.DP_crit
    P0plusDP = u[1]
    T = u[2]
    eps_g = u[3]
    V = u[4]
    rho_m = u[5]
    tot_m = u[8]
    tot_w = u[9]
    tot_c = u[10]

    P = P0plusDP + P_lit - param.P_lit_0

    m_h20 = tot_w / tot_m
    m_co2 = tot_c / tot_m

    eps_x = crystal_fraction_eps_x(composition, T, P, m_h20, m_co2)
    m_eq_max = exsolve_meq(composition, P, T, 0.0)

    # MT's new stuff
    eps_m0 = 1 - eps_x

    m_h2o_melt = tot_w / (V * rho_m * eps_m0)
    m_co2_melt = tot_c / (V * rho_m * eps_m0)

    C_co2_sat = exsolve3(composition, P, T, m_h2o_melt)[1]

    out[1] = eps_x
    out[2] = eps_x / (1 - eps_g) - 0.8
    out[3] = if sw.eruption == 0
        (P - P_lit) - DP_crit
    else
        -DP_crit
    end
    out[4] = if sw.eruption == 1
        P_lit - P
    else
        -DP_crit
    end
    out[5] = eps_x - 0.5
    out[6] = m_h2o_melt - m_eq_max
    out[7] = -(P0plusDP - param.P_lit_0 + DP_crit)
    return out[8] = m_co2_melt - C_co2_sat
end

"""
    affect!(int, idx, sw::SW{Int8}, param::Param{Float64}, param_saved_var::ParamSaved{Float64}, param_IC_Finder::ParamICFinder{Float64})

Re-initialize the condition when the event happens. This function modifies the current state of the integrator (`int.u`) when a particular event occurs during the simulation. The function adjusts various parameters based on the current state of the integrator and the custom parameters that were passed in.

# Arguments:
- `int`: The current state of the integrator. It's format is from the DifferentialEquations.jl package
- `idx`: The index of the event that caused the function to be called.
- `sw`: A custom parameter used to control simulation behavior.
- `param`: A custom parameter containing physical constants and other model parameters.
- `param_saved_var`: A custom parameter used to store values from the previous time step.
- `param_IC_Finder`: A custom parameter used to control the behavior of the IC_Finder function.

The arguments `int` and `idx` are from the DifferentialEquations.jl package. These argument formats are specific to the DifferentialEquations.jl package.
"""
function affect!(
    int,
    idx,
    sw::SW{Int8},
    param::Param{Float64},
    param_saved_var::ParamSaved{Float64},
    param_IC_Finder::ParamICFinder{Float64},
    erupt_saved::EruptSaved{Float64},
)
    println("*event idx: ", idx)
    composition = param.composition
    storeTime = param_saved_var.storeTime
    storeTemp = param_saved_var.storeTemp
    storeTemp = storeTemp[storeTime .< int.t]
    storeTime = storeTime[storeTime .< int.t]
    param_saved_var.storeTime = storeTime
    param_saved_var.storeTemp = storeTemp

    if param.dP_lit_dt_0 == 0
        temp_P_lit = 0.0
    else
        if int.t <= abs(param.P_lit_drop_max / param.dP_lit_dt_0)
            temp_P_lit = param.dP_lit_dt_0 * int.t
        else
            temp_P_lit = -param.P_lit_drop_max
        end
    end
    P_0 = int.u[1] + temp_P_lit

    m_h2o = int.u[9] / int.u[8]
    m_co2 = int.u[10] / int.u[8]
    eps_x0 = crystal_fraction_eps_x(composition, int.u[2], P_0, m_h2o, m_co2)

    if idx == 3 && eps_x0 < 0.5
        sw.eruption = 1
        rho_g0 = eos_g_rho_g(P_0, int.u[2])
        record_erupt_start(int.t, int.u[3], eps_x0, int.u[5], int.u[6], rho_g0, erupt_saved)
        println("reached critical pressure and need to start an eruption,  time: ", int.t)
    elseif idx == 4
        sw.eruption = 0
        record_erupt_end(int.t, erupt_saved, param)
        println("If it just finished an eruption...  time: ", int.t)
    elseif idx == 6 || idx == 8
        phase_here = param_saved_var.phase
        println(
            "starting ic finder for conversion of phase,  time: $(int.t), phase_here: $phase_here",
        )
        eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder(
            composition,
            int.u[9],
            int.u[10],
            int.u[8],
            P_0,
            int.u[2],
            int.u[4],
            int.u[5],
            param_IC_Finder,
        )

        param_saved_var.phase = phase
        if phase_here != phase
            println("1st try in IC Finder successful")
            int.u[3] = eps_g_temp
            int.u[7] = X_co2_temp
            C_co2 = C_co2_temp
        else
            println("trying new IC parameters...")
            param_IC_Finder.max_count = 150
            eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder(
                composition,
                int.u[9],
                int.u[10],
                int.u[8],
                P_0,
                int.u[2],
                int.u[4],
                int.u[5],
                param_IC_Finder,
            )
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
                param_IC_Finder.Tol = param_IC_Finder.Tol * 0.1
                eps_g_temp, X_co2_temp, C_co2_temp, phase = IC_Finder(
                    composition,
                    int.u[9],
                    int.u[10],
                    int.u[8],
                    P_0,
                    int.u[2],
                    int.u[4],
                    int.u[5],
                    param_IC_Finder,
                )
                param_saved_var.phase = phase
                ## change back to initial Tol
                param_IC_Finder.Tol = param_IC_Finder.Tol * 10
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
