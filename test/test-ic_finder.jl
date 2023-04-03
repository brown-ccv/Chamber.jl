include("ic_finder-silicic-data.jl")
@testset "IC_Finder-silicic" begin
    eps_g, X_co2, mco2_diss, phase = IC_Finder(
        composition, M_h2o, M_co2, M_tot, P, T, V, rho_m, mm_co2, mm_h2o, param_IC
    )
    @test (round(eps_g0; digits=17), X_co20, mco2_diss0, phase0) ==
        (round(eps_g; digits=17), X_co2, mco2_diss, phase)
end

include("ic_finder-mafic-data.jl")
@testset "IC_Finder-mafic" begin
    eps_g, X_co2, mco2_diss, phase = IC_Finder(
        composition, M_h2o, M_co2, M_tot, P, T, V, rho_m, mm_co2, mm_h2o, param_IC
    )
    @test (eps_g0, X_co20, mco2_diss0, phase0) == (eps_g, X_co2, mco2_diss, phase)
end
