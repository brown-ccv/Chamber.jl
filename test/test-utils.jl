@testset "constructer of sw" begin
    sw = SW()
    @test all([sw.heat_cond == 1, sw.visc_relax == 1, sw.eruption == 0])
end