function write_csv(df::DataFrame, path::String)::Nothing
    number_of_data = size(df, 1)
    println("number_of_data: $number_of_data")
    names = ["P+dP", "T", "eps_g", "V", "rho_m", "rho_x", "X_CO2", "total_mass", "total_mass_H2O", "total_mass_CO2"]
    rename!(df, ["value$i" => names[i] for i in 1:10])
    CSV.write("$path/out.csv", df)
    nothing
end