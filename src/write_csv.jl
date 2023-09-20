function extend_with_nan(arr::Vector{T}, n::Int64)::Vector{T} where {T<:Float64}
    return [arr; fill(NaN, n - length(arr))]
end

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
    df[!, "eps_x"] = crystal_fraction_eps_x.(repeat([composition], number_of_data), df[:, "T"], df[:, "P+dP"], m_h2o, m_co2)

    n = length(erupt_saved.time)
    df_erupt = DataFrame(
        time_of_eruption = erupt_saved.time,
        duration_of_eruption = extend_with_nan(erupt_saved.duration, n),
        mass_erupted = extend_with_nan(erupt_saved.mass, n),
        volume_erupted = extend_with_nan(erupt_saved.volume, n)
    )
    CSV.write("$path/out.csv", df)
    CSV.write("$path/eruptions.csv", df_erupt)
    return nothing
end
