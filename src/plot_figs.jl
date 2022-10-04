using DataFrames, CSV, Plots

function plot_figs(csv_file, path)
    ENV["GKSwstype"] = "100"  # magic environmental variable for Plots
    df  = CSV.File(csv_file) |> DataFrame
    p1 = plot(df[!,1], df[!,2], xaxis="Time", yaxis="P+dP", label="Julia", linewidth=2, marker=(:x,3,Plots.stroke(2)))
    p2 = plot(df[!,1], df[!,3], xaxis="Time", yaxis="T", label="Julia", linewidth=2, marker=(:x,3,Plots.stroke(2)))
    p3 = plot(df[!,1], df[!,4], xaxis="Time", yaxis="eps_g", label="Julia", linewidth=2, marker=(:x,3,Plots.stroke(2)))
    p4 = plot(df[!,1], df[!,5], xaxis="Time", yaxis="V", label="Julia", linewidth=2, marker=(:x,3,Plots.stroke(2)))
    p5 = plot(df[!,1], df[!,8], xaxis="Time", yaxis="X_CO2", label="Julia", linewidth=2, marker=(:x,3,Plots.stroke(2)))
    p6 = plot(df[!,1], df[!,9], xaxis="Time", yaxis="total_Mass", label="Julia", linewidth=2, marker=(:x,3,Plots.stroke(2)))

    savefig(p1, "$path/1_P.png")
    savefig(p2, "$path/2_T.png")
    savefig(p3, "$path/3_eps_g.png")
    savefig(p4, "$path/4_V.png")
    savefig(p5, "$path/5_X_CO2.png")
    savefig(p6, "$path/6_tot_Mass.png")
end
