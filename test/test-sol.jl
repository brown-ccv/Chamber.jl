using DelimitedFiles, DataFrames

solve_exp = Float64.(readdlm("out.csv", ',')[2:end, :])
df = chamber(Silicic(), 1e9, 0.2, 0.04, 0.001, -3.3, 8e3)

@testset "solve_ode" begin
    for i in 1:11
        println("col$i maximum diff: ", maximum!([1.0], abs.(solve_exp[:, i] .- Matrix(df)[:, i])))
        @test all(abs.(solve_exp[:, i] .- Matrix(df)[:, i]) .< abs(solve_exp[end, i]*1e-9))
    end
end
