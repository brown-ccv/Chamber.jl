function write_csv(df::DataFrame, erupt_saved::EruptSaved{Float64}, path::String, composition::Union{Silicic,Mafic})::Nothing
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
    rename!(df, ["value$i" => name for (i, name) in enumerate(names)])
    # add eps_x column
    m_h2o = df[:, "total_mass_H2O"] ./ df[:, "total_mass"]
    m_co2 = df[:, "total_mass_CO2"] ./ df[:, "total_mass"]
    crystal_fraction_eps_x′(T, P, mH2O, mCO2) = crystal_fraction_eps_x(composition, T, P, mH2O, mCO2)
    df[!, "eps_x"] = crystal_fraction_eps_x′.(df[:, "T"], df[:, "P+dP"], m_h2o, m_co2)

    n = length(erupt_saved.time)
    df_erupt = DataFrame(
        time_of_eruption = erupt_saved.time,
        duration_of_eruption = [erupt_saved.duration; fill(NaN, n - length(erupt_saved.duration))],
        mass_erupted = [erupt_saved.mass; fill(NaN, n - length(erupt_saved.mass))],
        volume_erupted = [erupt_saved.volume; fill(NaN, n - length(erupt_saved.volume))]
    )
    CSV.write("$path/out.csv", df)
    CSV.write("$path/eruptions.csv", df_erupt)
    return nothing
end
