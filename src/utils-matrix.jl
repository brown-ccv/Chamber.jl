"""
This is for a11 and a12
    a11: rho, drho_dP, V, dV_dP
    a12: rho, drho_dT, V, dV_dT
"""
a1x(rho, drho_dX, V, dV_dX) = (1/rho)*drho_dX + (1/V)*dV_dX
a13(rho, drho_deps_g) = (1/rho)*drho_deps_g

# a14 = 0

# conservation of water mass
"""
This is for a21 and a22
    a21: P
    a22: T
"""
a2x(eps_m, dm_eq_dX, m_eq, deps_x_dX, dV_dX, V, drho_m_dX, rho_m, X_co2, m_g, eps_g, mm_h2o, drho_g_dX, rho_g) = 
eps_m*dm_eq_dX-
m_eq*deps_x_dX+
m_eq*eps_m*dV_dX/V+
m_eq*eps_m*drho_m_dX/rho_m+
(1-X_co2)/m_g*eps_g*mm_h2o*drho_g_dX/rho_m+
(1-X_co2)/m_g*eps_g*mm_h2o*rho_g*dV_dX/(rho_m*V)

a23(m_eq, X_co2, m_g, rho_g, mm_h2o, rho_m) = -m_eq+(1-X_co2)/m_g*rho_g*mm_h2o/rho_m
a24(eps_m, dm_eq_dX_co2, m_g, eps_g, rho_g, mm_h2o, rho_m, X_co2, mm_co2) = 
real(eps_m*dm_eq_dX_co2-
1/m_g*eps_g*rho_g*mm_h2o/rho_m-
(1-X_co2)*eps_g*rho_g*mm_h2o*(mm_co2-mm_h2o)/(Complex(m_g)^2*rho_m))

# conservation of (total) enthalpy - all divided by rc*T*V
a31(drc_dP, rc, dV_dP, V, L_m, eps_x, drho_x_dP, T, rho_x, deps_x_dP, L_e, dm_eq_dP, rho_m, eps_m, m_eq, drho_m_dP) = 
    drc_dP/(rc)+
    dV_dP/V+
    L_m*(-eps_x*drho_x_dP/(rc*T)-
    rho_x*deps_x_dP/(rc*T)-
    rho_x*eps_x*dV_dP/(rc*V*T))+
    L_e*(-dm_eq_dP*rho_m*eps_m/(rc*T)-
    m_eq*eps_m*drho_m_dP/(rc*T)+
    m_eq*rho_m*deps_x_dP/(rc*T)-
    m_eq*rho_m*eps_m*dV_dP/(rc*V*T))

a32(drc_dT, rc, T, dV_dT, V, L_m, eps_x, drho_x_dT, rho_x, deps_x_dT, L_e, dm_eq_dT, rho_m, eps_m, m_eq, drho_m_dT) = 
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

a33(rho_g, c_g, rho_m, c_m, rc, L_e, m_eq, T) = (rho_g*c_g-rho_m*c_m)/rc +
    L_e*m_eq*rho_m/(rc*T)
a34(L_e, rho_m, eps_m, dm_eq_dX_co2, rc, T) = -L_e*rho_m*eps_m*dm_eq_dX_co2/(rc*T)

    # conservation of CO2 mass
"""
This is for a41 and a42
    a41: P
    a42: T
"""
a4x(eps_m, dC_co2_dX, C_co2, deps_x_dX, dV_dX, V, drho_m_dX, rho_m, X_co2, m_g, eps_g, mm_co2, drho_g_dX, rho_g) = 
    eps_m*dC_co2_dX-
    C_co2*deps_x_dX+
    C_co2*eps_m*dV_dX/V+
    C_co2*eps_m*drho_m_dX/rho_m+
    X_co2/m_g*eps_g*mm_co2*drho_g_dX/rho_m+
    X_co2/m_g*eps_g*rho_g*mm_co2*dV_dX/(rho_m*V)

a43(C_co2, X_co2, m_g, rho_g, mm_co2, rho_m) = -C_co2+
    X_co2/m_g*rho_g*mm_co2/rho_m

a44(eps_m, dC_co2_dX_co2, m_g, eps_g, rho_g, mm_co2, rho_m, X_co2, mm_h2o) = real(eps_m*dC_co2_dX_co2+
    1/m_g*eps_g*rho_g*mm_co2/rho_m-
    X_co2*eps_g*rho_g*mm_co2*(mm_co2-mm_h2o)/(Complex(m_g)^2*rho_m))



b1(Mdot_in, Mdot_out, rho, V, P_loss, rho_x, rho_m, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, a11, dP_lit_dt) = 
    (Mdot_in - Mdot_out)/(rho*V)-
    P_loss-
    (rho_x-rho_m)/rho*deps_x_dmh2o_t*dM_h2o_t_dt-
    (rho_x-rho_m)/rho*deps_x_dmco2_t*dM_co2_t_dt-
    a11*dP_lit_dt

    # conservation of water mass
b2(Mdot_v_in, Mdot_v_out, rho_m, V, m_eq, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, eps_m, P_loss, X_co2, m_g, rho_g, eps_g, mm_h2o, a21, dP_lit_dt) = 
    (Mdot_v_in - Mdot_v_out)/(rho_m*V)-
    m_eq*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)-
    m_eq*eps_m*P_loss-
    (1-X_co2)/m_g*rho_g/rho_m*eps_g*mm_h2o*P_loss-
    a21*dP_lit_dt
    # conservation of (total) enthalpy
b3(Hdot_in, Hdot_out, rc, T, V, rho_x, c_x, rho_m, c_m, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, L_m, L_e, m_eq, P_loss, eps_x, eps_m, a31, dP_lit_dt) = 
    (Hdot_in - Hdot_out)/(rc*T*V)-
    1/(rc*V*T)*((rho_x*c_x-rho_m*c_m)*T*V*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt))+
    L_m*rho_x*V/(rc*V*T)*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)+
    L_e*m_eq*rho_m*V/(rc*V*T)*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)+
    P_loss*(-1+L_m*rho_x*eps_x*V/(rc*V*T)+
    L_e*m_eq*rho_m*eps_m*V/(rc*V*T))-
    a31*dP_lit_dt
    # conservation of CO2 mass
b4(Mdot_c_in, Mdot_c_out, rho_m, V, C_co2, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, eps_m, P_loss, X_co2, m_g, rho_g, eps_g, mm_co2, a41, dP_lit_dt) = 
    (Mdot_c_in - Mdot_c_out)/(rho_m*V)-
    C_co2*(deps_x_dmh2o_t*dM_h2o_t_dt+deps_x_dmco2_t*dM_co2_t_dt)-
    C_co2*eps_m*P_loss-
    (X_co2)/m_g*rho_g/rho_m*eps_g*mm_co2*P_loss-
    a41*dP_lit_dt

    # if phase == 3
    #     # set up matrices to solve using Cramer"s rule
    #     A          = [a11 a12 a13 a14; a21 a22 a23 a24; a31 a32 a33 a34; a41 a42 a43 a44]
    #     A_P        = [b1  a12 a13 a14; b2  a22 a23 a24; b3  a32 a33 a34; b4  a42 a43 a44]
    #     A_T        = [a11 b1  a13 a14; a21 b2  a23 a24; a31 b3  a33 a34; a41 b4  a43 a44]
    #     A_eps_g    = [a11 a12 b1  a14; a21 a22 b2  a24; a31 a32 b3  a34; a41 a42 b4  a44]
    #     A_X_co2    = [a11 a12 a13 b1 ; a21 a22 a23 b2 ; a31 a32 a33 b3 ; a41 a42 a43 b4 ]

    #     det_A          = det(A)
    #     dDP_dt         = det(A_P)/det_A
    #     dT_dt          = det(A_T)/det_A
    #     deps_g_dt      = det(A_eps_g)/det_A
    #     dX_co2_dt      = det(A_X_co2)/det_A
    #     dP_dt          = dDP_dt+param["dP_lit_dt"]
    #     dV_dt          = dV_dP*dP_dt + dV_dT*dT_dt + V*P_loss
    #     drho_m_dt      = drho_m_dP*dP_dt + drho_m_dT*dT_dt
    #     drho_x_dt      = drho_x_dP*dP_dt + drho_x_dT*dT_dt
    # else

    #     A          = [a11 a12; a31 a32]
    #     A_P        = [b1  a12; b3  a32]
    #     A_T        = [a11 b1 ; a31 b3]
        
    #     det_A          = det(A)
    #     dDP_dt         = det(A_P)/det_A
    #     dT_dt          = det(A_T)/det_A
    #     deps_g_dt      = 0
            
    #     dP_dt          = dDP_dt+param["dP_lit_dt"]
    #     dV_dt          = dV_dP*dP_dt + dV_dT*dT_dt + V*P_loss
    #     drho_m_dt      = drho_m_dP*dP_dt + drho_m_dT*dT_dt
    #     drho_x_dt      = drho_x_dP*dP_dt + drho_x_dT*dT_dt
    #     dX_co2_dt      = 0
    # end

function build_matrix(phase::Int, rho::Float64, drho_dP::Float64, V::Float64, dV_dP::Float64, drho_dT::Float64, dV_dT::Float64, drc_dP::Float64, rc::Float64, L_m::Float64, eps_x::Float64, drho_x_dP::Float64, T::Float64, 
        rho_x::Float64, deps_x_dP::Float64, L_e::Float64, dm_eq_dP::Float64, rho_m::Float64, eps_m::Float64, m_eq::Float64, drho_m_dP::Float64, drc_dT::Float64, drho_x_dT::Float64, deps_x_dT::Float64, dm_eq_dT::Float64, 
        drho_m_dT::Float64, Mdot_in::Float64, Mdot_out::Float64, P_loss::Float64, deps_x_dmh2o_t::Float64, dM_h2o_t_dt::Float64, deps_x_dmco2_t::Float64, dM_co2_t_dt::Float64, dP_lit_dt::Float64, 
        Hdot_in::Float64, Hdot_out::Float64, c_x::Float64, c_m::Float64, drho_deps_g::Float64, X_co2::Float64, m_g::Float64, eps_g::Float64, mm_h2o::Float64, drho_g_dP::Float64, rho_g::Float64, drho_g_dT::Float64, dm_eq_dX_co2::Float64, 
        mm_co2::Float64, c_g::Float64, dC_co2_dP::Float64, C_co2::Float64, dC_co2_dT::Float64, dC_co2_dX_co2::Float64, Mdot_v_in::Float64, Mdot_v_out::Float64, Mdot_c_in::Float64, Mdot_c_out::Float64)

    a110 = a1x(rho, drho_dP, V, dV_dP)
    a120 = a1x(rho, drho_dT, V, dV_dT)
    a310 = a31(drc_dP, rc, dV_dP, V, L_m, eps_x, drho_x_dP, T, rho_x, deps_x_dP, L_e, dm_eq_dP, rho_m, eps_m, m_eq, drho_m_dP)
    a320 = a32(drc_dT, rc, T, dV_dT, V, L_m, eps_x, drho_x_dT, rho_x, deps_x_dT, L_e, dm_eq_dT, rho_m, eps_m, m_eq, drho_m_dT)
    b10 = b1(Mdot_in, Mdot_out, rho, V, P_loss, rho_x, rho_m, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, a110, dP_lit_dt)
    b30 = b3(Hdot_in, Hdot_out, rc, T, V, rho_x, c_x, rho_m, c_m, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, L_m, L_e, m_eq, P_loss, eps_x, eps_m, a310, dP_lit_dt)
    if phase != 3
        A = [a110 a120; a310 a320]
        b = [b10, b30]
        return [A, b]
    elseif phase == 3
        a130 = a13(rho, drho_deps_g)
        a140 = 0
        a210 = a2x(eps_m, dm_eq_dP, m_eq, deps_x_dP, dV_dP, V, drho_m_dP, rho_m, X_co2, m_g, eps_g, mm_h2o, drho_g_dP, rho_g)
        a220 = a2x(eps_m, dm_eq_dT, m_eq, deps_x_dT, dV_dT, V, drho_m_dT, rho_m, X_co2, m_g, eps_g, mm_h2o, drho_g_dT, rho_g)
        a230 = a23(m_eq, X_co2, m_g, rho_g, mm_h2o, rho_m)
        a240 = a24(eps_m, dm_eq_dX_co2, m_g, eps_g, rho_g, mm_h2o, rho_m, X_co2, mm_co2)
        a330 = a33(rho_g, c_g, rho_m, c_m, rc, L_e, m_eq, T)
        a340 = a34(L_e, rho_m, eps_m, dm_eq_dX_co2, rc, T)
        a410 = a4x(eps_m, dC_co2_dP, C_co2, deps_x_dP, dV_dP, V, drho_m_dP, rho_m, X_co2, m_g, eps_g, mm_co2, drho_g_dP, rho_g)
        a420 = a4x(eps_m, dC_co2_dT, C_co2, deps_x_dT, dV_dT, V, drho_m_dT, rho_m, X_co2, m_g, eps_g, mm_co2, drho_g_dT, rho_g)
        a430 = a43(C_co2, X_co2, m_g, rho_g, mm_co2, rho_m)
        a440 = a44(eps_m, dC_co2_dX_co2, m_g, eps_g, rho_g, mm_co2, rho_m, X_co2, mm_h2o)

        b20 = b2(Mdot_v_in, Mdot_v_out, rho_m, V, m_eq, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, eps_m, P_loss, X_co2, m_g, rho_g, eps_g, mm_h2o, a210, dP_lit_dt)
        b40 = b4(Mdot_c_in, Mdot_c_out, rho_m, V, C_co2, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, eps_m, P_loss, X_co2, m_g, rho_g, eps_g, mm_co2, a410, dP_lit_dt)

        A = [a110 a120 a130 a140; a210 a220 a230 a240; a310 a320 a330 a340; a410 a420 a430 a440]
        b = [b10, b20, b30, b40]
        return [A, b]
    else
        error("phase not found")
    end
end
