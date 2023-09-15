function plot_sol(df, x_axis_col, y_axis_col)
    return plot(
        df[:, x_axis_col],
        df[:, y_axis_col];
        xaxis=x_axis_col,
        yaxis=y_axis_col,
        linewidth=2,
        marker=(:x, 3, Plots.stroke(2)),
    )
end

function plot_sol_year_unit(df, xlabel, ylabel)
    # Units of outputs: Time -> sec, Presure -> Pa, Volume -> m³
    divisors = Dict(
        "time" => 31536000,
        "P+dP" => 1e6,
        "V" => 1e9,
        "total_mass" => 1e12,
    )
    labels = Dict(
        "time" => "Time (yr)",
        "P+dP" => "Pressure (MPa)",
        "T" => "Temperature (K)",
        "V" => "Volume (km³)",
        "total_mass" => "Total mass (×10¹² kg)",
    )
    convert_unit(label) = df[:, label] / get(divisors, label, 1.0)
    getaxis(label) = get(labels, label, label)

    return plot(
        convert_unit(xlabel), convert_unit(ylabel);
        xaxis=getaxis(xlabel),
        yaxis=getaxis(ylabel),
        legend=false,
        linewidth=2,
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

function plot_combined_fig(df::DataFrame, path::String, title::LaTeXString)::Nothing
    ENV["GKSwstype"] = "100"  # magic environmental variable for Plots
    plot_names = ["P+dP", "T", "eps_g", "V", "X_CO2", "total_mass"]
    plts = [plot_sol_year_unit(df, "time", l) for l in plot_names]
    combined_plot = plot(plts..., layout=(2, 3), size=(900, 500), dpi=200, margin = 3Plots.mm, plot_title=title, plot_titlefontsize=13)
    savefig(combined_plot, joinpath(path, "combined_plots.png"))
    return nothing
end

function plot_ϵx(storeTime::Vector{Float64}, storeEps_x::Vector{Float64}, path::String)::Nothing
    ENV["GKSwstype"] = "100"  # magic environmental variable for Plots
    savefig(plot(
            storeTime, storeEps_x;
            xaxis="time",
            yaxis="eps_x",
            label="Chamber.jl",
            linewidth=2,
            marker=(:x, 3, Plots.stroke(2))),
        "$path/eps_x.png")
    return nothing
end

function plot_dual_axis(df::DataFrame, storeTime::Vector{Float64}, storeEps_x::Vector{Float64}, path::String)::Nothing
    color1 = RGB(20/255, 125/255, 204/255)
    color2 = RGB(255/255, 100/255, 0/255)
    ENV["GKSwstype"] = "100"  # magic environmental variable for Plots

    common_opt1 = (
        xaxis="Time (yr)",
        linewidth=3,
        legend=:false,
        guidefontsize=18,
        tickfontsize=12,
        linecolor=color1,
        fg_color_guide=color1,
        fg_color_text=color1,
    )

    common_opt2 = (
        linewidth=3,
        legend=:false,
        guidefontsize=18,
        tickfontsize=12,
        linecolor=color2,
        fg_color_guide=color2,
        fg_color_text=color2,
    )

    plot(
        df[!, 1]/31536000,
        df[!, 2]/1e6;
        common_opt1...,
        yaxis="Pressure (MPa)",
        )
    p1 = plot!(
        twinx(),
        df[!, 1]/31536000,
        df[!, 3];
        common_opt2...,
        yaxis="Temperature (K)",
    )

    plot(
        df[!, 1]/31536000,
        df[!, 4];
        common_opt1...,
        yaxis=L"𝜺_g",
        )
    p2 = plot!(
        twinx(),
        storeTime/31536000,
        storeEps_x;
        common_opt2...,
        yaxis=L"𝜺_X",
    )

    savefig(p1, joinpath(path, "P_T.png"))
    savefig(p2, joinpath(path, "eps_x_eps_g.png"))
    return nothing
end
