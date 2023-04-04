include("test-initial-data.jl")
@testset "find_liq" begin
    @test find_liq(Silicic(), h2o, co2, P, ini_eps_x) == Tl_silicic
    @test find_liq(Mafic(), h2o, co2, P, ini_eps_x) == Tl_mafic
end

@testset "crystal_fraction-silicic" begin
    @test (
        eps_x_s,
        deps_x_dP_s,
        deps_x_dT_s,
        deps_x_deps_g_s,
        deps_x_dmco2_t_s,
        deps_x_dmh2o_t_s,
    ) == values(crystal_fraction(Silicic(), T_s, P_s, mH2O_s, mCO2_s))
    @test crystal_fraction_eps_x(Silicic(), T_s, P_s, mH2O_s, mCO2_s) == eps_x_s
end

@testset "crystal_fraction-mafic" begin
    @test (
        eps_x_m,
        deps_x_dP_m,
        deps_x_dT_m,
        deps_x_deps_g_m,
        deps_x_dmco2_t_m,
        deps_x_dmh2o_t_m,
    ) == values(crystal_fraction(Mafic(), T_m, P_m, mH2O_m, mCO2_m))
    @test crystal_fraction_eps_x(Mafic(), T_m, P_m, mH2O_m, mCO2_m) == eps_x_m
end

@testset "parameter_melting_curve-silicic" begin
    @test (a_s, dadx_s, dady_s, dadz_s, b_s, dbdx_s, dbdy_s, dbdz_s, c_s, dcdx_s, dcdy_s, dcdz_s) == values(parameters_melting_curve(Silicic(), h2o, co2, P))
end

@testset "parameter_melting_curve-mafic" begin
    @test (a_m, dadx_m, dady_m, dadz_m, b_m, dbdx_m, dbdy_m, dbdz_m) == values(parameters_melting_curve(Mafic(), h2o, co2, P))
end
