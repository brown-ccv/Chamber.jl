include("exsolve-silicic-data.jl")
@testset "exsolve-silicic" begin
    @test (; meq, dmeqdP, dmeqdT, dmeqdXco2, C_co2, dC_co2dP, dC_co2dT, dC_co2dXco2) ==
        exsolve(composition, P, T, X_co2)
    @test meq == exsolve_meq(composition, P, T, X_co2)
end

include("exsolve-mafic-data.jl")
@testset "exsolve-mafic" begin
    @test (; meq, dmeqdP, dmeqdT, dmeqdXco2, C_co2, dC_co2dP, dC_co2dT, dC_co2dXco2) ==
        exsolve(composition, P, T, X_co2)
    @test meq == exsolve_meq(composition, P, T, X_co2)
end

include("exsolve3-silicic-data.jl")
@testset "exsolve3-silicic" begin
    @test (; C_co2, X_co2) == exsolve3(composition, P, T, m_h2o_melt)
end

include("exsolve3-mafic-data.jl")
@testset "exsolve3-mafic" begin
    @test (; C_co2, X_co2) == exsolve3(composition, P, T, m_h2o_melt)
end
