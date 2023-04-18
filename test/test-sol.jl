using DelimitedFiles, DataFrames

solve_exp = Float64.(readdlm("out.csv", ',')[2:end, :])
df = chamber(Silicic(), 1e9, 0.2, 0.04, 0.001, -3.3, 8e3)

@testset "solve_ode" begin
    @test all(abs.(solve_exp .- Matrix(df)) .< 0.00005)
end
