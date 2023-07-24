function plot_sol(df, x_axis_col, y_axis_col)
    return plot(
        df[:, x_axis_col],
        df[:, y_axis_col];
        xaxis=x_axis_col,
        yaxis=y_axis_col,
        label="Chamber.jl",
        linewidth=2,
        marker=(:x, 3, Plots.stroke(2)),
    )
end

function plot_figs(df::DataFrame, path::String)::Nothing
    ENV["GKSwstype"] = "100"  # magic environmental variable for Plots
    savefig(plot_sol(df, "time", "P+dP"), "$path/1_P.png")
    savefig(plot_sol(df, "time", "T"), "$path/2_T.png")
    savefig(plot_sol(df, "time", "eps_g"), "$path/3_eps_g.png")
    savefig(plot_sol(df, "time", "V"), "$path/4_V.png")
    savefig(plot_sol(df, "time", "X_CO2"), "$path/5_X_CO2.png")
    savefig(plot_sol(df, "time", "total_mass"), "$path/6_tot_Mass.png")
    return nothing
end
