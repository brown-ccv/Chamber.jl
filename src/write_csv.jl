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

function write_ϵx_csv(storeTime::Vector{Float64}, storeEps_x::Vector{Float64}, path::String)::Nothing
    df = DataFrame(Time = storeTime, eps_x = storeEps_x)
    CSV.write("$path/eps_x.csv", df)
    return nothing
end
