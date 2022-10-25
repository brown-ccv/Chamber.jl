
s = """a11 = a1x(rho, drho_dP, V, dV_dP)
a12 = a1x(rho, drho_dT, V, dV_dT)
a31 = a31(drc_dP, rc, dV_dP, V, L_m, eps_x, drho_x_dP, T, rho_x, deps_x_dP, L_e, dm_eq_dP, rho_m, eps_m, m_eq, drho_m_dP)
a32 = a32(drc_dT, rc, T, dV_dT, V, L_m, eps_x, drho_x_dT, rho_x, deps_x_dT, L_e, dm_eq_dT, rho_m, eps_m, m_eq, drho_m_dT)
b1 = b1(Mdot_in, Mdot_out, rho, V, P_loss, rho_x, rho_m, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, a11, dP_lit_dt)
b3 = b3(Hdot_in, Hdot_out, rc, T, V, rho_x, c_x, rho_m, c_m, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, L_m, L_e, m_eq, P_loss, eps_x, eps_m, a31, dP_lit_dt)
a13 = a13(rho, drho_deps_g)
a21 = a2x(eps_m, dm_eq_dP, m_eq, deps_x_dP, dV_dP, V, drho_m_dP, rho_m, X_co2, m_g, eps_g, mm_h2o, drho_g_dP, rho_g)
a22 = a2x(eps_m, dm_eq_dT, m_eq, deps_x_dT, dV_dT, V, drho_m_dT, rho_m, X_co2, m_g, eps_g, mm_h2o, drho_g_dT, rho_g)
a23 = a23(m_eq, X_co2, m_g, rho_g, mm_h2o, rho_m)
a24 = a24(eps_m, dm_eq_dX_co2, m_g, eps_g, rho_g, mm_h2o, rho_m, X_co2, mm_co2)
a33 = a33(rho_g, c_g, rho_m, c_m, rc, L_e, m_eq, T)
a34 = a34(L_e, rho_m, eps_m, dm_eq_dX_co2, rc, T)
a41 = a4x(eps_m, dC_co2_dP, C_co2, deps_x_dP, dV_dP, V, drho_m_dP, rho_m, X_co2, m_g, eps_g, mm_co2, drho_g_dP, rho_g)
a42 = a4x(eps_m, dC_co2_dT, C_co2, deps_x_dT, dV_dT, V, drho_m_dT, rho_m, X_co2, m_g, eps_g, mm_co2, drho_g_dT, rho_g)
a43 = a43(C_co2, X_co2, m_g, rho_g, mm_co2, rho_m)
a44 = a44(eps_m, dC_co2_dX_co2, m_g, eps_g, rho_g, mm_co2, rho_m, X_co2, mm_h2o)
b2 = b2(Mdot_v_in, Mdot_v_out, rho_m, V, m_eq, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, eps_m, P_loss, X_co2, m_g, rho_g, eps_g, mm_h2o, a21, dP_lit_dt)
b4 = b4(Mdot_c_in, Mdot_c_out, rho_m, V, C_co2, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, eps_m, P_loss, X_co2, m_g, rho_g, eps_g, mm_co2, a41, dP_lit_dt)
"""

s1 = split(s, "\n")
s2 = [split(s, "(") for s in s1][1:end-1]
s3 = [s[2][1:end-1] for s in s2]
s4 = []
for s in s3
    for t in split(s, ", ")
        push!(s4, t)
    end
end
s5 = join(unique(s4), ", ")