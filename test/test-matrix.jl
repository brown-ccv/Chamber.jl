include("test-matrix-a-data.jl")
@testset "Matrix A" begin
    @test a11 == a1x_f(rho, drho_dP, V, dV_dP)
    @test a12 == a1x_f(rho, drho_dT, V, dV_dT)
    @test a13 == a13_f(rho, drho_deps_g)

    @test a21 == a21_f(
        eps_m,
        dm_eq_dP,
        C_h2o,
        deps_x_dP,
        dV_dP,
        V,
        drho_m_dP,
        rho_m,
        X_co2,
        m_g,
        eps_g,
        mm_h2o,
        drho_g_dP,
        rho_g,
    )
    @test a22 == a22_f(
        eps_m,
        dm_eq_dT,
        C_h2o,
        deps_x_dT,
        dV_dT,
        V,
        drho_m_dT,
        rho_m,
        X_co2,
        m_g,
        eps_g,
        mm_h2o,
        drho_g_dT,
        rho_g,
    )
    @test a23 == a23_f(C_h2o, X_co2, m_g, rho_g, mm_h2o, rho_m)
    @test a24 == a24_f(eps_m, dm_eq_dX_co2, m_g, eps_g, rho_g, mm_h2o, rho_m, X_co2, mm_co2)

    @test a31 == a31_f(
        drc_dP,
        rc,
        dV_dP,
        V,
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
        C_h2o,
        drho_m_dP,
    )
    @test a32 == a32_f(
        drc_dT,
        rc,
        T,
        dV_dT,
        V,
        L_m,
        eps_x,
        drho_x_dT,
        rho_x,
        deps_x_dT,
        L_e,
        dm_eq_dT,
        rho_m,
        eps_m,
        C_h2o,
        drho_m_dT,
    )
    @test a33 == a33_f(rho_g, c_g, rho_m, c_m, rc, L_e, C_h2o, T)
    @test a34 == a34_f(L_e, rho_m, eps_m, dm_eq_dX_co2, rc, T)

    @test a41 == a41_f(
        eps_m,
        dC_co2_dP,
        C_co2,
        deps_x_dP,
        dV_dP,
        V,
        drho_m_dP,
        rho_m,
        X_co2,
        m_g,
        eps_g,
        mm_co2,
        drho_g_dP,
        rho_g,
    )
    @test a42 == a42_f(
        eps_m,
        dC_co2_dT,
        C_co2,
        deps_x_dT,
        dV_dT,
        V,
        drho_m_dT,
        rho_m,
        X_co2,
        m_g,
        eps_g,
        mm_co2,
        drho_g_dT,
        rho_g,
    )
    @test a43 == a43_f(C_co2, X_co2, m_g, rho_g, mm_co2, rho_m)
    @test a44 ==
        a44_f(eps_m, dC_co2_dX_co2, m_g, eps_g, rho_g, mm_co2, rho_m, X_co2, mm_h2o)
end

include("test-matrix-b-data.jl")
@testset "Matrix B" begin
    @test b1 == b1_f(
        Mdot_in,
        Mdot_out,
        rho,
        V,
        P_loss,
        rho_x,
        rho_m,
        deps_x_dmh2o_t,
        dM_h2o_t_dt,
        deps_x_dmco2_t,
        dM_co2_t_dt,
    )
    @test b2 == b2_f(
        Mdot_v_in,
        Mdot_v_out,
        rho_m,
        V,
        C_h2o,
        deps_x_dmh2o_t,
        dM_h2o_t_dt,
        deps_x_dmco2_t,
        dM_co2_t_dt,
        eps_m,
        P_loss,
        X_co2,
        m_g,
        rho_g,
        eps_g,
        mm_h2o,
    )
    @test b3 == b3_f(
        Hdot_in,
        Hdot_out,
        rc,
        T,
        V,
        rho_x,
        c_x,
        rho_m,
        c_m,
        deps_x_dmh2o_t,
        dM_h2o_t_dt,
        deps_x_dmco2_t,
        dM_co2_t_dt,
        L_m,
        L_e,
        C_h2o,
        P_loss,
        eps_x,
        eps_m,
    )
    @test b4 == b4_f(
        Mdot_c_in,
        Mdot_c_out,
        rho_m,
        V,
        C_co2,
        deps_x_dmh2o_t,
        dM_h2o_t_dt,
        deps_x_dmco2_t,
        dM_co2_t_dt,
        eps_m,
        P_loss,
        X_co2,
        m_g,
        rho_g,
        eps_g,
        mm_co2,
    )
end

phase = 2.0
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
    C_h2o,
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

@testset "Build Matrix phase=2" begin
    @test A == [a11 a12; a31 a32]
    @test b == [b1, b3]
end

phase = 3.0
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
    C_h2o,
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
x = A \ b
@testset "Build Matrix phase=3" begin
    @test A == [a11 a12 a13 a14; a21 a22 a23 a24; a31 a32 a33 a34; a41 a42 a43 a44]
    @test b == [b1, b2, b3, b4]
    @test x[1] == dP_dt
    @test round(x[2]; digits=24) == round(dT_dt; digits=24)
    @test round(x[3]; digits=26) == round(deps_g_dt; digits=26)
    @test round(x[4]; digits=25) == round(dX_co2_dt; digits=25)
    @test dV_dP * dP_dt + dV_dT * dT_dt + V * P_loss == dV_dt
    @test drho_m_dP * dP_dt + drho_m_dT * dT_dt == drho_m_dt
    @test drho_x_dP * dP_dt + drho_x_dT * dT_dt == drho_x_dt
end
