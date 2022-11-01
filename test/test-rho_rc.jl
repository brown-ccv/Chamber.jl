@testset "rho-rc" begin
    @test rho == rho_f(eps_m=eps_m, eps_g=eps_g, eps_x=eps_x, rho_m=rho_m, rho_g=rho_g, rho_x=rho_x)
    @test drho_dP == drho_dX_f(eps_m=eps_m, eps_g=eps_g, eps_x=eps_x, drho_m_dX=drho_m_dP, drho_g_dX=drho_g_dP, drho_x_dX=drho_x_dP, rho_x=rho_x, rho_m=rho_m, deps_x_dX=deps_x_dP) 
    @test drho_dT == drho_dX_f(eps_m=eps_m, eps_g=eps_g, eps_x=eps_x, drho_m_dX=drho_m_dT, drho_g_dX=drho_g_dT, drho_x_dX=drho_x_dT, rho_x=rho_x, rho_m=rho_m, deps_x_dX=deps_x_dT)
    @test rc == rc_f(rho_x=rho_x, eps_x=eps_x, c_x=c_x, rho_m=rho_m, eps_m=eps_m, c_m=c_m, rho_g=rho_g, eps_g=eps_g, c_g=c_g) 

    @test drc_dP == drc_dX_f(eps_x=eps_x, c_x=c_x, drho_x_dX=drho_x_dP, eps_g=eps_g, c_g=c_g, drho_g_dX=drho_g_dP, eps_m=eps_m, c_m=c_m, drho_m_dX=drho_m_dP, rho_x=rho_x, rho_m=rho_m, deps_x_dX=deps_x_dP)
    @test drc_dT == drc_dX_f(eps_x=eps_x, c_x=c_x, drho_x_dX=drho_x_dT, eps_g=eps_g, c_g=c_g, drho_g_dX=drho_g_dT, eps_m=eps_m, c_m=c_m, drho_m_dX=drho_m_dT, rho_x=rho_x, rho_m=rho_m, deps_x_dX=deps_x_dT)
    @test [rho, drho_dP, drho_dT, drho_deps_g, rc, drc_dP, drc_dT] == build_rho_rc(eps_m, eps_g, eps_x, rho_m, rho_g, rho_x, drho_m_dP, drho_g_dP, drho_x_dP, drho_m_dT, drho_g_dT, drho_x_dT, c_x, c_m, c_g, deps_x_dP, deps_x_dT)
end











