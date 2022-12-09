@testset "constructer of sw" begin
    sw = SW()
    @test all([sw.heat_cond == 1, sw.visc_relax == 1, sw.eruption == 0])
end

@testset "build_mdot_in" begin
    @test build_mdot_in(false, 2400.0, -3.3, 2.1582e8, 1079.8869633986928) == 38.14210301577417
end