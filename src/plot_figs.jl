function plot_sol(df, x_axis_col, y_axis_col)
    t = df[:, x_axis_col]
    y = df[:, y_axis_col]
    return plot(
        t,
        y;
        xaxis=x_axis_col,
        yaxis=y_axis_col,
        label="Chamber.jl",
        linewidth=2,
        marker=(:x, 3, Plots.stroke(2)),
    )
end

function plot_figs(df::DataFrame, path::String)::Nothing
    ENV["GKSwstype"] = "100"  # magic environmental variable for Plots
    p1 = plot_sol(df, "time", "P+dP")
    p2 = plot_sol(df, "time", "T")
    p3 = plot_sol(df, "time", "eps_g")
    p4 = plot_sol(df, "time", "V")
    p5 = plot_sol(df, "time", "X_CO2")
    p6 = plot_sol(df, "time", "total_mass")

    savefig(p1, "$path/1_P.png")
    savefig(p2, "$path/2_T.png")
    savefig(p3, "$path/3_eps_g.png")
    savefig(p4, "$path/4_V.png")
    savefig(p5, "$path/5_X_CO2.png")
    savefig(p6, "$path/6_tot_Mass.png")
    return nothing
end
