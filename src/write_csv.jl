function write_csv(df::DataFrame, erupt_saved::EruptSaved{Float64}, path::String)::Nothing
    number_of_data = size(df, 1)
    println("number_of_data: $number_of_data")
    names = [
        "P+dP",
        "T",
        "eps_g",
        "V",
        "rho_m",
        "rho_x",
        "X_CO2",
        "total_mass",
        "total_mass_H2O",
        "total_mass_CO2",
    ]
    rename!(df, ["timestamp" => "time"])
    rename!(df, ["value$i" => names[i] for i in 1:10])
    df_erupt = DataFrame(
        reduce(
            hcat,
            [erupt_saved.time, erupt_saved.duration, erupt_saved.mass, erupt_saved.volume],
        ),
        ["time of eruption", "duration of eruption", "mass erupted", "volume erupted"],
    )
    CSV.write("$path/out.csv", df)
    CSV.write("$path/eruptions.csv", df_erupt)
    return nothing
end
