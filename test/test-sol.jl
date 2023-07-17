using DelimitedFiles, DataFrames

solve_exp = Float64.(readdlm("out.csv", ',')[2:end, :])
df = chamber(Silicic(), 1e9, 0.2, 0.04, 0.001, -3.3, 8e3)

println("1: ", solve_exp[:, 1] .- Matrix(df)[:, 1])
println("2: ", solve_exp[:, 2] .- Matrix(df)[:, 2])
println("3: ", solve_exp[:, 3] .- Matrix(df)[:, 3])
println("4: ", solve_exp[:, 4] .- Matrix(df)[:, 4])
println("5: ", solve_exp[:, 5] .- Matrix(df)[:, 5])
println("6: ", solve_exp[:, 6] .- Matrix(df)[:, 6])
println("7: ", solve_exp[:, 7] .- Matrix(df)[:, 7])
println("8: ", solve_exp[:, 8] .- Matrix(df)[:, 8])
println("9: ", solve_exp[:, 9] .- Matrix(df)[:, 9])
println("10: ", solve_exp[:, 10] .- Matrix(df)[:, 10])
@testset "solve_ode" begin
    @test all(abs.(solve_exp .- Matrix(df)) .< 1)
    @test all(abs.(solve_exp[:, 1] .- Matrix(df)[:, 1]) .< 100.0)
    @test all(abs.(solve_exp[:, 2] .- Matrix(df)[:, 2]) .< abs(solve_exp[end, 2]*1e-5))
    @test all(abs.(solve_exp[:, 3] .- Matrix(df)[:, 3]) .< abs(solve_exp[end, 3]*1e-5))
    @test all(abs.(solve_exp[:, 4] .- Matrix(df)[:, 4]) .< abs(solve_exp[end, 4]*1e-5))
    @test all(abs.(solve_exp[:, 5] .- Matrix(df)[:, 5]) .< abs(solve_exp[end, 5]*1e-5))
    @test all(abs.(solve_exp[:, 6] .- Matrix(df)[:, 6]) .< abs(solve_exp[end, 6]*1e-5))
    @test all(abs.(solve_exp[:, 7] .- Matrix(df)[:, 7]) .< abs(solve_exp[end, 7]*1e-5))
    @test all(abs.(solve_exp[:, 8] .- Matrix(df)[:, 8]) .< abs(solve_exp[end, 8]*1e-5))
    @test all(abs.(solve_exp[:, 9] .- Matrix(df)[:, 9]) .< abs(solve_exp[end, 9]*1e-5))
    @test all(abs.(solve_exp[:, 10] .- Matrix(df)[:, 10]) .< abs(solve_exp[end, 10]*1e-5))
end
