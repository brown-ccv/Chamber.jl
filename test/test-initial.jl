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
