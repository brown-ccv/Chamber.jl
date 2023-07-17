using DelimitedFiles, DataFrames

solve_exp = Float64.(readdlm("out.csv", ',')[2:end, :])
df = chamber(Silicic(), 1e9, 0.2, 0.04, 0.001, -3.3, 8e3)

@testset "solve_ode" begin
    @test all(abs.(solve_exp .- Matrix(df)) .< 0.00005)
end

solve_exp1 = Float64.(readdlm("out_co2_1e-3.csv", ',')[2:end, :])
solve_exp2 = Float64.(readdlm("out_co2_2e-3.csv", ',')[2:end, :])
output_dirname = chamber(Mafic(), 1e9, 0.2, 0.01, [0.001, 0.002], -3.3, 8e3)
dir1 = joinpath(pwd(), output_dirname, "vol0.2_h2o0.01_gas0.001_vfr-3.3_dep8000.0")
dir2 = joinpath(pwd(), output_dirname, "vol0.2_h2o0.01_gas0.002_vfr-3.3_dep8000.0")
solve_out1 = Float64.(readdlm(joinpath(dir1, "out.csv"), ',')[2:end, :])
solve_out2 = Float64.(readdlm(joinpath(dir2, "out.csv"), ',')[2:end, :])
@testset "solve_ode array" begin
    @test isdir(dir1)
    @test isdir(dir2)
    @test isfile(joinpath(dir1, "out.csv"))
    @test isfile(joinpath(dir2, "out.csv"))
    for i in 1:10
        @test all(abs.(solve_exp1[:, i] .- solve_out1[:, i]) .< solve_exp1[end, i]*1e-3)
        @test all(abs.(solve_exp2[:, i] .- solve_out2[:, i]) .< solve_exp2[end, i]*1e-3)
    end
end
