include("ic_finder-silicic-data.jl")
@testset "IC_Finder-silicic" begin
    # Case 1
    eps_g, X_co2, mco2_diss, phase = IC_Finder(
        composition, M_h2o, M_co2, M_tot, P, T, V, rho_m, param_IC
    )
    @test (round(eps_g0; digits=17), X_co20, mco2_diss0, phase0) ==
        (round(eps_g; digits=17), X_co2, mco2_diss, phase)

    # Case 2
    eps_g, X_co2, mco2_diss, phase = IC_Finder(
        composition, M_h2o_2, M_co2_2, M_tot_2, P_2, T_2, V_2, rho_m_2, param_IC
    )
    @test (eps_g0_2, X_co20_2, mco2_diss0_2, phase0_2) == (eps_g, X_co2, mco2_diss, phase)

    # Case 3
    param_IC.max_count = 150
    eps_g, X_co2, mco2_diss, phase = IC_Finder(
        composition, M_h2o_3, M_co2_3, M_tot_3, P_3, T_3, V_3, rho_m_3, param_IC
    )
    @test (eps_g0_3, X_co20_3, mco2_diss0_3, phase0_3) == (eps_g, X_co2, mco2_diss, phase)
    param_IC.max_count = 100
end

include("ic_finder-mafic-data.jl")
@testset "IC_Finder-mafic" begin
    # Case 1
    eps_g, X_co2, mco2_diss, phase = IC_Finder(
        composition, M_h2o, M_co2, M_tot, P, T, V, rho_m, param_IC
    )
    @test (eps_g0, X_co20, mco2_diss0, phase0) == (eps_g, X_co2, mco2_diss, phase)

    # Case 2
    eps_g, X_co2, mco2_diss, phase = IC_Finder(
        composition, M_h2o_2, M_co2_2, M_tot_2, P_2, T_2, V_2, rho_m_2, param_IC
    )
    @test (eps_g0_2, X_co20_2, mco2_diss0_2, phase0_2) == (eps_g, X_co2, mco2_diss, phase)

    # Case 3
    eps_g, X_co2, mco2_diss, phase = IC_Finder(
        composition, M_h2o_3, M_co2_3, M_tot_3, P_3, T_3, V_3, rho_m_3, param_IC
    )
    @test (eps_g0_3, X_co20_3, mco2_diss0_3, phase0_3) == (eps_g, X_co2, mco2_diss, phase)

    # Case 4
    param_IC.max_count = 150
    eps_g, X_co2, mco2_diss, phase = IC_Finder(
        composition, M_h2o_4, M_co2_4, M_tot_4, P_4, T_4, V_4, rho_m_4, param_IC
    )
    @test (eps_g0_4, X_co20_4, mco2_diss0_4, phase0_4) == (eps_g, X_co2, mco2_diss, phase)
    param_IC.max_count = 100
end
