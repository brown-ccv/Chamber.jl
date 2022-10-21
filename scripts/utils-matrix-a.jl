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
a21(eps_m, dm_eq_dX, m_eq, deps_x_dX, dV_dX, V, drho_m_dX, rho_m, X_co2, m_g, eps_g, mm_h2o, drho_g_dX, rho_g) = 
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
