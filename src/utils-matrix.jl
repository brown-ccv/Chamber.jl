# Matrix A
# a1x: conservation of (total) mass
"""
This is for a11 and a12
    a11: rho, drho_dP, V, dV_dP
    a12: rho, drho_dT, V, dV_dT
"""
a1x_f(rho, drho_dX, V, dV_dX) = (1/rho)*drho_dX + (1/V)*dV_dX
a13_f(rho, drho_deps_g) = (1/rho)*drho_deps_g

# a2x: conservation of water mass
a21_f(eps_m, dm_eq_dP, m_eq, deps_x_dP, dV_dP, V, drho_m_dP, rho_m, X_co2, m_g, eps_g, mm_h2o, drho_g_dP, rho_g) =
eps_m*dm_eq_dP-
m_eq*deps_x_dP+
m_eq*eps_m*dV_dP/V+
m_eq*eps_m*drho_m_dP/rho_m+
(1-X_co2)/m_g*eps_g*mm_h2o*drho_g_dP/rho_m+
(1-X_co2)/m_g*eps_g*mm_h2o*rho_g*dV_dP/(rho_m*V)

a22_f(eps_m, dm_eq_dT, m_eq, deps_x_dT, dV_dT, V, drho_m_dT, rho_m, X_co2, m_g, eps_g, mm_h2o, drho_g_dT, rho_g) =
eps_m*dm_eq_dT-
m_eq*deps_x_dT-
m_eq*eps_m*dV_dT/V+
m_eq*eps_m*drho_m_dT/rho_m+
(1-X_co2)/m_g*eps_g*mm_h2o*drho_g_dT/rho_m+
(1-X_co2)/m_g*eps_g*mm_h2o*rho_g*dV_dT/(rho_m*V)

a23_f(m_eq, X_co2, m_g, rho_g, mm_h2o, rho_m) = -m_eq+(1-X_co2)/m_g*rho_g*mm_h2o/rho_m
a24_f(eps_m, dm_eq_dX_co2, m_g, eps_g, rho_g, mm_h2o, rho_m, X_co2, mm_co2) = 
real(eps_m*dm_eq_dX_co2-
1/m_g*eps_g*rho_g*mm_h2o/rho_m-
(1-X_co2)*eps_g*rho_g*mm_h2o*(mm_co2-mm_h2o)/(Complex(m_g)^2*rho_m))

# a3x: conservation of (total) enthalpy - all divided by rc*T*V
a31_f(drc_dP, rc, dV_dP, V, L_m, eps_x, drho_x_dP, T, rho_x, deps_x_dP, L_e, dm_eq_dP, rho_m, eps_m, m_eq, drho_m_dP) = 
    drc_dP/(rc)+
    dV_dP/V+
    L_m*(-eps_x*drho_x_dP/(rc*T)-
    rho_x*deps_x_dP/(rc*T)-
    rho_x*eps_x*dV_dP/(rc*V*T))+
    L_e*(-dm_eq_dP*rho_m*eps_m/(rc*T)-
    m_eq*eps_m*drho_m_dP/(rc*T)+
    m_eq*rho_m*deps_x_dP/(rc*T)-
    m_eq*rho_m*eps_m*dV_dP/(rc*V*T))

a32_f(drc_dT, rc, T, dV_dT, V, L_m, eps_x, drho_x_dT, rho_x, deps_x_dT, L_e, dm_eq_dT, rho_m, eps_m, m_eq, drho_m_dT) = 
    drc_dT/(rc)+
    1/T+
    dV_dT/V+
    L_m*(-eps_x*drho_x_dT/(rc*T)-
    rho_x*deps_x_dT/(rc*T)-
    rho_x*eps_x*dV_dT/(rc*T*V))+
    L_e*(-dm_eq_dT*rho_m*eps_m/(rc*T)-
    m_eq*eps_m*drho_m_dT/(rc*T)+
    m_eq*rho_m*deps_x_dT/(rc*T)-
    m_eq*rho_m*eps_m*dV_dT/(rc*T*V))

a33_f(rho_g, c_g, rho_m, c_m, rc, L_e, m_eq, T) = (rho_g*c_g-rho_m*c_m)/rc +
    L_e*m_eq*rho_m/(rc*T)
a34_f(L_e, rho_m, eps_m, dm_eq_dX_co2, rc, T) = -L_e*rho_m*eps_m*dm_eq_dX_co2/(rc*T)

# a4x: conservation of CO2 mass
a41_f(eps_m, dC_co2_dP, C_co2, deps_x_dP, dV_dP, V, drho_m_dP, rho_m, X_co2, m_g, eps_g, mm_co2, drho_g_dP, rho_g) = 
    eps_m*dC_co2_dP-
    C_co2*deps_x_dP+
    C_co2*eps_m*dV_dP/V+
    C_co2*eps_m*drho_m_dP/rho_m+
    X_co2/m_g*eps_g*mm_co2*drho_g_dP/rho_m+
    X_co2/m_g*eps_g*rho_g*mm_co2*dV_dP/(rho_m*V)

a42_f(eps_m, dC_co2_dT, C_co2, deps_x_dT, dV_dT, V, drho_m_dT, rho_m, X_co2, m_g, eps_g, mm_co2, drho_g_dT, rho_g) = 
    eps_m*dC_co2_dT-
    C_co2*deps_x_dT-
    C_co2*eps_m*dV_dT/V+
    C_co2*eps_m*drho_m_dT/rho_m+
    X_co2/m_g*eps_g*mm_co2*drho_g_dT/rho_m+
    X_co2/m_g*eps_g*rho_g*mm_co2*dV_dT/(rho_m*V)

a43_f(C_co2, X_co2, m_g, rho_g, mm_co2, rho_m) = -C_co2+
    X_co2/m_g*rho_g*mm_co2/rho_m

a44_f(eps_m, dC_co2_dX_co2, m_g, eps_g, rho_g, mm_co2, rho_m, X_co2, mm_h2o) = real(eps_m*dC_co2_dX_co2+
    1/m_g*eps_g*rho_g*mm_co2/rho_m-
    X_co2*eps_g*rho_g*mm_co2*(mm_co2-mm_h2o)/(Complex(m_g)^2*rho_m))


# Vector B
"""
    dM_X_t_dt_f(rho, V, Mdot_v_in, Mdot_v_out, m_h2o, Mdot_in, Mdot_out)

For `dM_h2o_t_dt` & `dM_co2_t_dt`

- Mdot_X_in: `Mdot_v_in` or `Mdot_c_in`
- Mdot_X_out: `Mdot_v_out` or `Mdot_c_out`
- m_X: `m_h2o` or `m_co2`
"""
dM_X_t_dt_f(rho, V, Mdot_X_in, Mdot_X_out, m_X, Mdot_in, Mdot_out) = 
    1/(rho*V)*((Mdot_X_in-Mdot_X_out)-m_X*(Mdot_in-Mdot_out))

# b1: conservation of (total) mass
b1_f(Mdot_in, Mdot_out, rho, V, P_loss, rho_x, rho_m, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, a11, dP_lit_dt) = 
    (Mdot_in - Mdot_out)/(rho*V)-
    P_loss-
    (rho_x-rho_m)/rho*deps_x_dmh2o_t*dM_h2o_t_dt-
    (rho_x-rho_m)/rho*deps_x_dmco2_t*dM_co2_t_dt-
    a11*dP_lit_dt

# b2: conservation of water mass
b2_f(Mdot_v_in, Mdot_v_out, rho_m, V, m_eq, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, eps_m, P_loss, X_co2, m_g, rho_g, eps_g, mm_h2o, a21, dP_lit_dt) = 
    (Mdot_v_in - Mdot_v_out)/(rho_m*V)-
    m_eq*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)-
    m_eq*eps_m*P_loss-
    (1-X_co2)/m_g*rho_g/rho_m*eps_g*mm_h2o*P_loss-
    a21*dP_lit_dt

# b3: conservation of (total) enthalpy
b3_f(Hdot_in, Hdot_out, rc, T, V, rho_x, c_x, rho_m, c_m, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, L_m, L_e, m_eq, P_loss, eps_x, eps_m, a31, dP_lit_dt) = 
    (Hdot_in - Hdot_out)/(rc*T*V)-
    1/(rc*V*T)*((rho_x*c_x-rho_m*c_m)*T*V*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt))+
    L_m*rho_x*V/(rc*V*T)*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)+
    L_e*m_eq*rho_m*V/(rc*V*T)*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)+
    P_loss*(-1+L_m*rho_x*eps_x*V/(rc*V*T)+
    L_e*m_eq*rho_m*eps_m*V/(rc*V*T))-
    a31*dP_lit_dt

# b4: conservation of CO2 mass
b4_f(Mdot_c_in, Mdot_c_out, rho_m, V, C_co2, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, eps_m, P_loss, X_co2, m_g, rho_g, eps_g, mm_co2, a41, dP_lit_dt) = 
    (Mdot_c_in - Mdot_c_out)/(rho_m*V)-
    C_co2*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)-
    C_co2*eps_m*P_loss-
    (X_co2)/m_g*rho_g/rho_m*eps_g*mm_co2*P_loss-
    a41*dP_lit_dt

function build_matrix(phase::Int64, rho::Float64, drho_dP::Float64, V::Float64, dV_dP::Float64, drho_dT::Float64, dV_dT::Float64, drc_dP::Float64, rc::Float64, L_m::Float64, eps_x::Float64, drho_x_dP::Float64, T::Float64, 
        rho_x::Float64, deps_x_dP::Float64, L_e::Float64, dm_eq_dP::Float64, rho_m::Float64, eps_m::Float64, m_eq::Float64, drho_m_dP::Float64, drc_dT::Float64, drho_x_dT::Float64, deps_x_dT::Float64, dm_eq_dT::Float64, 
        drho_m_dT::Float64, Mdot_in::Float64, Mdot_out::Float64, P_loss::Float64, deps_x_dmh2o_t::Float64, m_h2o::Float64, m_co2::Float64, deps_x_dmco2_t::Float64, dP_lit_dt::Float64, 
        Hdot_in::Float64, Hdot_out::Float64, c_x::Float64, c_m::Float64, drho_deps_g::Float64, X_co2::Float64, m_g::Float64, eps_g::Float64, mm_h2o::Float64, drho_g_dP::Float64, rho_g::Float64, drho_g_dT::Float64, dm_eq_dX_co2::Float64, 
        mm_co2::Float64, c_g::Float64, dC_co2_dP::Float64, C_co2::Float64, dC_co2_dT::Float64, dC_co2_dX_co2::Float64, Mdot_v_in::Float64, Mdot_v_out::Float64, Mdot_c_in::Float64, Mdot_c_out::Float64)

    a11 = a1x_f(rho, drho_dP, V, dV_dP)
    a12 = a1x_f(rho, drho_dT, V, dV_dT)
    a31 = a31_f(drc_dP, rc, dV_dP, V, L_m, eps_x, drho_x_dP, T, rho_x, deps_x_dP, L_e, dm_eq_dP, rho_m, eps_m, m_eq, drho_m_dP)
    a32 = a32_f(drc_dT, rc, T, dV_dT, V, L_m, eps_x, drho_x_dT, rho_x, deps_x_dT, L_e, dm_eq_dT, rho_m, eps_m, m_eq, drho_m_dT)
    dM_h2o_t_dt = dM_X_t_dt_f(rho, V, Mdot_v_in, Mdot_v_out, m_h2o, Mdot_in, Mdot_out)
    dM_co2_t_dt = dM_X_t_dt_f(rho, V, Mdot_c_in, Mdot_c_out, m_co2, Mdot_in, Mdot_out)
    b1 = b1_f(Mdot_in, Mdot_out, rho, V, P_loss, rho_x, rho_m, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, a11, dP_lit_dt)
    b3 = b3_f(Hdot_in, Hdot_out, rc, T, V, rho_x, c_x, rho_m, c_m, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, L_m, L_e, m_eq, P_loss, eps_x, eps_m, a31, dP_lit_dt)
    if phase != 3
        A = [a11 a12; a31 a32]
        b = [b1, b3]
        return [A, b]
    elseif phase == 3
        a13 = a13_f(rho, drho_deps_g)
        a14 = 0
        a21 = a21_f(eps_m, dm_eq_dP, m_eq, deps_x_dP, dV_dP, V, drho_m_dP, rho_m, X_co2, m_g, eps_g, mm_h2o, drho_g_dP, rho_g)
        a22 = a22_f(eps_m, dm_eq_dT, m_eq, deps_x_dT, dV_dT, V, drho_m_dT, rho_m, X_co2, m_g, eps_g, mm_h2o, drho_g_dT, rho_g)
        a23 = a23_f(m_eq, X_co2, m_g, rho_g, mm_h2o, rho_m)
        a24 = a24_f(eps_m, dm_eq_dX_co2, m_g, eps_g, rho_g, mm_h2o, rho_m, X_co2, mm_co2)
        a33 = a33_f(rho_g, c_g, rho_m, c_m, rc, L_e, m_eq, T)
        a34 = a34_f(L_e, rho_m, eps_m, dm_eq_dX_co2, rc, T)
        a41 = a41_f(eps_m, dC_co2_dP, C_co2, deps_x_dP, dV_dP, V, drho_m_dP, rho_m, X_co2, m_g, eps_g, mm_co2, drho_g_dP, rho_g)
        a42 = a42_f(eps_m, dC_co2_dT, C_co2, deps_x_dT, dV_dT, V, drho_m_dT, rho_m, X_co2, m_g, eps_g, mm_co2, drho_g_dT, rho_g)
        a43 = a43_f(C_co2, X_co2, m_g, rho_g, mm_co2, rho_m)
        a44 = a44_f(eps_m, dC_co2_dX_co2, m_g, eps_g, rho_g, mm_co2, rho_m, X_co2, mm_h2o)

        b2 = b2_f(Mdot_v_in, Mdot_v_out, rho_m, V, m_eq, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, eps_m, P_loss, X_co2, m_g, rho_g, eps_g, mm_h2o, a21, dP_lit_dt)
        b4 = b4_f(Mdot_c_in, Mdot_c_out, rho_m, V, C_co2, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, eps_m, P_loss, X_co2, m_g, rho_g, eps_g, mm_co2, a41, dP_lit_dt)

        A = [a11 a12 a13 a14; a21 a22 a23 a24; a31 a32 a33 a34; a41 a42 a43 a44]
        b = [b1, b2, b3, b4]
        return [A, b]
    else
        error("phase not found")
    end
end
