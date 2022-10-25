include("test-matrix-a-data.jl")
@testset "Matrix A" begin
    @test a110 == a1x(rho, drho_dP, V, dV_dP)
    @test a120 == a1x(rho, drho_dT, V, dV_dT)
    @test a130 == a13(rho, drho_deps_g)

    @test a210 == a2x(eps_m, dm_eq_dP, m_eq, deps_x_dP, dV_dP, V, drho_m_dP, rho_m, X_co2, m_g, eps_g, mm_h2o, drho_g_dP, rho_g)
    @test round(a220, digits=5) == round(a2x(eps_m, dm_eq_dT, m_eq, deps_x_dT, dV_dT, V, drho_m_dT, rho_m, X_co2, m_g, eps_g, mm_h2o, drho_g_dT, rho_g), digits=5)
    @test a230 == a23(m_eq, X_co2, m_g, rho_g, mm_h2o, rho_m)
    @test a240 == a24(eps_m, dm_eq_dX_co2, m_g, eps_g, rho_g, mm_h2o, rho_m, X_co2, mm_co2)

    @test a310 == a31(drc_dP, rc, dV_dP, V, L_m, eps_x, drho_x_dP, T, rho_x, deps_x_dP, L_e, dm_eq_dP, rho_m, eps_m, m_eq, drho_m_dP)
    @test a320 == a32(drc_dT, rc, T, dV_dT, V, L_m, eps_x, drho_x_dT, rho_x, deps_x_dT, L_e, dm_eq_dT, rho_m, eps_m, m_eq, drho_m_dT)
    @test a330 == a33(rho_g, c_g, rho_m, c_m, rc, L_e, m_eq, T)
    @test a340 == a34(L_e, rho_m, eps_m, dm_eq_dX_co2, rc, T)

    @test a410 == a4x(eps_m, dC_co2_dP, C_co2, deps_x_dP, dV_dP, V, drho_m_dP, rho_m, X_co2, m_g, eps_g, mm_co2, drho_g_dP, rho_g)
    @test round(a420, digits=6) == round(a4x(eps_m, dC_co2_dT, C_co2, deps_x_dT, dV_dT, V, drho_m_dT, rho_m, X_co2, m_g, eps_g, mm_co2, drho_g_dT, rho_g), digits=6)
    @test a430 == a43(C_co2, X_co2, m_g, rho_g, mm_co2, rho_m)
    @test a440 == a44(eps_m, dC_co2_dX_co2, m_g, eps_g, rho_g, mm_co2, rho_m, X_co2, mm_h2o)
end

include("test-matrix-b-data.jl")
@testset "Matrix B" begin
    @test b10 == b1(Mdot_in, Mdot_out, rho, V, P_loss, rho_x, rho_m, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, a110, dP_lit_dt)
    @test b20 == b2(Mdot_v_in, Mdot_v_out, rho_m, V, m_eq, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, eps_m, P_loss, X_co2, m_g, rho_g, eps_g, mm_h2o, a210, dP_lit_dt)
    @test b30 == b3(Hdot_in, Hdot_out, rc, T, V, rho_x, c_x, rho_m, c_m, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, L_m, L_e, m_eq, P_loss, eps_x, eps_m, a310, dP_lit_dt)
    @test b40 == b4(Mdot_c_in, Mdot_c_out, rho_m, V, C_co2, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, eps_m, P_loss, X_co2, m_g, rho_g, eps_g, mm_co2, a410, dP_lit_dt)
end

phase = 2
A, b = build_matrix(phase, rho, drho_dP, V, dV_dP, drho_dT, dV_dT, drc_dP, rc, L_m, eps_x, drho_x_dP, T, 
        rho_x, deps_x_dP, L_e, dm_eq_dP, rho_m, eps_m, m_eq, drho_m_dP, drc_dT, drho_x_dT, deps_x_dT, dm_eq_dT, 
        drho_m_dT, Mdot_in, Mdot_out, P_loss, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, dP_lit_dt, 
        Hdot_in, Hdot_out, c_x, c_m, drho_deps_g, X_co2, m_g, eps_g, mm_h2o, drho_g_dP, rho_g, drho_g_dT, dm_eq_dX_co2, 
        mm_co2, c_g, dC_co2_dP, C_co2, dC_co2_dT, dC_co2_dX_co2, Mdot_v_in, Mdot_v_out, Mdot_c_in, Mdot_c_out)

@testset "Build Matrix phase=2" begin
    @test A == [a110 a120; a310 a320]
    @test b == [b10, b30]
end

phase = 3
A, b = build_matrix(phase, rho, drho_dP, V, dV_dP, drho_dT, dV_dT, drc_dP, rc, L_m, eps_x, drho_x_dP, T, 
        rho_x, deps_x_dP, L_e, dm_eq_dP, rho_m, eps_m, m_eq, drho_m_dP, drc_dT, drho_x_dT, deps_x_dT, dm_eq_dT, 
        drho_m_dT, Mdot_in, Mdot_out, P_loss, deps_x_dmh2o_t, dM_h2o_t_dt, deps_x_dmco2_t, dM_co2_t_dt, dP_lit_dt, 
        Hdot_in, Hdot_out, c_x, c_m, drho_deps_g, X_co2, m_g, eps_g, mm_h2o, drho_g_dP, rho_g, drho_g_dT, dm_eq_dX_co2, 
        mm_co2, c_g, dC_co2_dP, C_co2, dC_co2_dT, dC_co2_dX_co2, Mdot_v_in, Mdot_v_out, Mdot_c_in, Mdot_c_out)

@testset "Build Matrix phase=3" begin
    @test A == [a110 a120 a130 a140; a210 a220 a230 a240; a310 a320 a330 a340; a410 a420 a430 a440]
    @test b == [b10, b20, b30, b40]
end