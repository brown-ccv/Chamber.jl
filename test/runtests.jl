using Chamber
using Test
include("matrix-a-data.jl")

@testset "GLQ_points_weights test" begin
    @test Chamber.GLQ_points_weights_hard(5) == [0 0]
    A = [1 1.0000000000000000 -0.5773502691896257; 2 1.0000000000000000 0.5773502691896257]
    @test Chamber.GLQ_points_weights_hard(2) == [A[:,3], A[:,2]]
    A = [1 0.3607615730481386 0.6612093864662645;2 0.3607615730481386 -0.6612093864662645;3 0.4679139345726910 -0.2386191860831969;4 0.4679139345726910 0.2386191860831969;5 0.1713244923791704 -0.9324695142031521;6 0.1713244923791704 0.9324695142031521]
    @test Chamber.GLQ_points_weights_hard(6) == [A[:,3], A[:,2]]
end

@testset "gas_heat_capacity test" begin
    @test Chamber.gas_heat_capacity(0) == 0
    @test Chamber.gas_heat_capacity(0.5) == 1978.5523133967436
    @test Chamber.gas_heat_capacity(1) == 1200.0
end

@testset "initial test" begin
    @test Chamber.eos_g(2.1582e8, 1000.0).rho_g == 502.5694116183761
    @test Chamber.eos_g(2.1582e8, 1000.0).drho_g_dP == 1.3610427024362434e-6
    @test Chamber.eos_g(2.1582e8, 1000.0).drho_g_dT == -0.649635523875751
    @test round(Chamber.eos_g(0.0, 1000.0).drho_g_dT, digits=5) == 4.79215
end

# @testset "export functions" begin
#     @test !isnothing(make_param)
#     @test !isnothing(GLQ_points_weights_hard)
# end

@testset "export functions Chamber" begin
    @test !isnothing(Chamber.make_param)
    @test !isnothing(Chamber.GLQ_points_weights_hard)
end

include("test-utils.jl")
include("test-rho_rc.jl")
