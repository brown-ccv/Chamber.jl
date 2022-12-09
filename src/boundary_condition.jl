function p_loss(visc_relax::Int8, P::Float64, P_lit::Float64, eta_r::Float64)::Float64
    if !(visc_relax in [0, 1])
        error("visc_relax must be 0 or 1.")
    else
        if visc_relax == 1
            P_loss = (P - P_lit)/eta_r 
        elseif visc_relax == 0
            P_loss = 0.0
        end
        return P_loss
    end
end

"""
    boundary_conditions_new(P::Float64, T::Float64, V::Float64, rho_m::Float64, rho_x::Float64, c::Float64, sw::SW{Int8}, T_in::Float64, M_h2o::Float64, M_co2::Float64, total_Mass::Float64, param::Param{Float64}, param_saved_var::ParamSaved{Float64})

# Arguments
`P`: Pressure (Pa)
`T`: Temperature (K)
`V`: chamber volume (m^3)
`rho_m`: density of melt
`rho_x`: density of crystal of magma
`c`: heat of magma
`sw`: eruption/cooling_module/viscous_relaxation control
`T_in`: Temperature
`M_h2o`: total mess of H2O in the magma
`M_co2`: total mess of CO2 in the magma
`total_Mass`: total mess of magma chamber
"""
function boundary_conditions_new(P::Float64, T::Float64, V::Float64, rho_m::Float64, rho_x::Float64, c::Float64, sw::SW{Int8}, T_in::Float64, M_h2o::Float64, M_co2::Float64, total_Mass::Float64, param::Param{Float64}, param_saved_var::ParamSaved{Float64})::Vector{Float64}
    P_lit = param.P_lit
    tot_h2o_frac_in = param.tot_h2o_frac_in
    tot_co2_frac_in = param.tot_co2_frac_in
    Mdot_in_pass = param.Mdot_in_pass
    Mdot_out_pass = param.Mdot_out_pass
    DP_crit = param.DP_crit
    Q_out_old = param.Q_out_old

    eps_g_in       = 0.0
    X_co2_in       = 0.0
    if param.fluxing
        X_co2_in = param.XCO2_in
    end
    rho_g_in = eos_g_rho_g(P, T_in)

    eps_x_in = crystal_fraction_eps_x(param.composition,T_in,P_lit,tot_h2o_frac_in,tot_co2_frac_in)

    rho_in = (1-eps_g_in-eps_x_in)*rho_m + eps_g_in*rho_g_in + eps_x_in*rho_x
    c_g_in = gas_heat_capacity(X_co2_in)
    c_in   = ((1-eps_g_in-eps_x_in)*rho_m*param.c_m + eps_g_in*rho_g_in*c_g_in + eps_x_in*rho_x*param.c_x)/rho_in
    if param.fluxing
        c_in = c_g_in
    end
    Mdot_in        = Mdot_in_pass
    Mdot_v_in      = tot_h2o_frac_in*Mdot_in
    Hdot_in        = c_in*T_in*Mdot_in
    Mdot_c_in      = tot_co2_frac_in*Mdot_in

    # set outflow conditions
    if sw.eruption == 0
        Mdot_out   = 0.0
        Mdot_v_out = 0.0
        Mdot_c_out = 0.0
    elseif sw.eruption == 1
        Mdot_out   = Mdot_out_pass
        Mdot_v_out = M_h2o/total_Mass*Mdot_out_pass
        Mdot_c_out = M_co2/total_Mass*Mdot_out_pass
    else
        println("eruption not specified")
    end

    a                    = real(Complex(V/(4*pi/3))^(1/3))   # chamber radius (m)
    cc                   = 10*a   # outer shell radius (m)
    dr                   = 0.1*a
    if ~isreal(a)
        println("here")
        a                    = real((V/(4*pi/3))^(1/3))   # chamber radius (m)
        cc                   = 10*a   # outer shell radius (m)
        dr                   = 0.1*a
    end

    if sw.heat_cond == 1
        # heat loss
        Q_out = heat_conduction_chamberCH(param.maxn,a,cc,dr,param.kappa,param.rho_r,param.c_r,param.Tb,param_saved_var)

    elseif sw.heat_cond == 0
        Q_out = 0.0
    else
        println("heat_cond not specified")
    end

    if isnan(Q_out)
        Q_out=Q_out_old
    else
        Q_out_old=Q_out
    end

    Hdot_out       = c*T*Mdot_out + Q_out

    # viscous relaxation
    quadpts,weights = GLQ_points_weights_hard(param.GLQ_n)
    b     = a + cc
    quadpts_r = (b-a)/2*quadpts .+ (a+b)/2

    Trt = heat_conduction_chamber_profileCH(param.maxn,a,cc,quadpts_r,param.kappa,param.Tb,param_saved_var)
    if param.rheol == "new"
        A = param.A  # material-dependent constant for viscosity law (Pa s)
        B = param.B  # molar gas constant (J/mol/K)
        G = param.G  # activation energy for creep (J/mol)
        eta_rt     = A*exp(G/B/Trt)
    elseif param.rheol == "old"
        nn = param.nn
        AA = param.AA
        G = param.G
        M = param.M
        dev_stress = DP_crit
        eta_rt     = (dev_stress^(1-nn)/AA)*exp(G/M/Trt)
    end
    I          = (b-a)/2*sum(weights*(eta_rt/(quadpts_r).^4))
    eta_r      = 3*a^3*I
    P_loss = p_loss(sw.visc_relax, P, P_lit, eta_r)

    return [Mdot_in, Mdot_out, Mdot_v_in, Mdot_v_out, Mdot_c_in, Mdot_c_out, Hdot_in, Hdot_out, P_loss, eta_r]
end