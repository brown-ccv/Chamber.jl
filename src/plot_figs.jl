function plot_sol(df, x_axis_col, y_axis_col)
    return plot(
        df[:, x_axis_col],
        df[:, y_axis_col];
        xaxis=x_axis_col,
        yaxis=y_axis_col,
        linewidth=2,
        legend=false,
        marker=(:x, 3, Plots.stroke(2)),
    )
end

function plot_figs(df::DataFrame, path::String)::Nothing
    ENV["GKSwstype"] = "100"  # magic environmental variable for Plots
    savefig′(label, path) = savefig(plot_sol(df, "time", label), path)
    labels = ["P+dP", "T", "eps_g", "V", "X_CO2", "total_mass"]
    paths = ["$path/$i-$l.png" for (i, l) in enumerate(labels)]
    broadcast(savefig′, labels, paths)
    return nothing
end
