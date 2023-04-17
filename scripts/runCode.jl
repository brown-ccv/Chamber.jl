using Chamber
include("./solver_methods.jl")

"""
    chamber(composition::Union{Silicic,Mafic}, end_time::Float64, log_volume_km3::Float64, InitialConc_H2O::Float64, InitialConc_CO2::Float64, log_vfr::Float64, depth::Float64, methods::Dict=methods, method::String="Tsit5", odesetting=OdeSetting(), ini_eps_x::Float64=0.15, rheol::String="old")

Simulate the eruption of a volcano using a model for the frequency of eruptions of upper crustal magma chambers based on Degruyter and Huber (2014).

# Arguments
- `composition`: Either a `Silicic` or a `Mafic` object that specifies the composition of the magma.
- `end_time`: Simulation period in seconds.
- `log_volume_km3`: The estimated volume of the volcano chamber in logarithmic scale.
- `InitialConc_H2O`: The initial water content of the magma, as a fraction of the total mass.
- `InitialConc_CO2`: The initial carbon dioxide content of the magma, as a fraction of the total mass.
- `log_vfr`: The logarithm of the volume flux rate, in cubic meters per second.
- `depth`: Estimate of depth of volcano chamber in meters.

# Returns
- A `DataFrame` containing the solution with columns: timestamp, P+dP, T, eps_g, V, rho_m, rho_x, X_CO2, total_mass, total_mass_H2O, and total_mass_CO2.

# Output
- `out.csv`: A CSV file containing the solution with headers: timestamp, P+dP, T, eps_g, V, rho_m, rho_x, X_CO2, total_mass, total_mass_H2O, and total_mass_CO2.
- Figures of the following variables plotted against time: P+dP, T, eps_g, V, X_CO2, total_mass.

# References
- W. Degruyter and C. Huber. A model for eruption frequency of upper crustal silicic magma chambers. Earth Planet. Sci. Lett. (2014).

# Examples
```
# Run a simulation with silicic magma chamber
julia> composition = Silicic()
julia> end_time = 3e9
julia> log_volume_km3 = 0.2
julia> InitialConc_H2O = 0.04
julia> InitialConc_CO2 = 0.001
julia> log_vfr = -3.3
julia> depth = 8e3

julia> chamber(composition, end_time, log_volume_km3, InitialConc_H2O, InitialConc_CO2, log_vfr, depth)

# Run a simulation with mafic magma chamber
julia> composition = Mafic()
julia> end_time = 3e9
julia> log_volume_km3 = 0.2
julia> InitialConc_H2O = 0.01
julia> InitialConc_CO2 = 0.001
julia> log_vfr = -3.3
julia> depth = 8e3

julia> chamber(composition, end_time, log_volume_km3, InitialConc_H2O, InitialConc_CO2, log_vfr, depth)
```

"""
function chamber(
    composition::Union{Silicic,Mafic},
    end_time::Float64,
    log_volume_km3::Float64,
    InitialConc_H2O::Float64,
    InitialConc_CO2::Float64,
    log_vfr::Float64,
    depth::Float64,
    methods::Dict=methods,
    method::String="Tsit5",
    odesetting=OdeSetting(),
    ini_eps_x::Float64=0.15,
    rheol::String="old",
)
    datetime = get_timestamp()
    composition_str = string(typeof(composition))
    path = joinpath(pwd(), "$(datetime)_$composition_str")
    mkdir(path)
    io = open("$path/$datetime.log", "w+")
    logger = SimpleLogger(io)
    global_logger(logger)

    to = get_timer("share")
    @timeit to "chamber" begin
        rc = rheol_composition_dict[composition_str]
        param = Param{Float64}(;
            composition=composition,
            rheol=rheol,
            rho_m0=rc.rho_m0,
            rho_x0=rc.rho_x0,
            c_m=rc.c_m,
            c_x=rc.c_x,
            L_m=rc.L_m,
        )
        r = rheol_dict[rheol]
        if rheol == "new"
            param.A, param.B, param.G = r.A, r.B, r.G
        elseif rheol == "old"
            param.nn, param.AA, param.G, param.M = r.nn, r.AA, r.G, r.M
        end

        c = ConstantValues()
        param_IC_Finder = ParamICFinder()
        param_saved_var = ParamSaved()
        sw = SW()

        # Initial temperature and viscosity profile around chamber
        param_saved_var.storeSumk = zeros(param.maxn)
        param_saved_var.storeSumk_2 = zeros(param.maxn)
        param_saved_var.storeSumk_old = zeros(param.maxn)
        param_saved_var.storeSumk_2_old = zeros(param.maxn)

        volume_km3 = 10^log_volume_km3                             # range of volume in km3
        range_radius = 1000 * (volume_km3 / (4 * pi / 3))^(1 / 3)  # range of radius in m
        V_0 = 4 * pi / 3 * range_radius^3                          # initial volume of the chamber (m^3)

        # thermal gradient
        param.Tb = c.T_surface + c.T_gradient * depth        # background temperature crust (K)

        # lithostatic pressure
        P_0 = param.rho_r * c.grav_acc * depth   # initial chamber pressure (Pa)
        param.P_lit = P_0
        param.P_lit_0 = P_0
        if param.single_eruption
            P_0 = P_0 + param.DP_crit
        end

        T_0 = find_liq(composition, InitialConc_H2O, InitialConc_CO2, P_0, ini_eps_x)

        T_in = T_0 + 50        # Temperature of inflowing magma (K)
        param.T_in = T_in

        # set the mass inflow rate
        param.Mdot_in_pass = build_mdot_in(param.fluxing, rc.rho_m0, log_vfr, P_0, T_in)

        rho_g0 = eos_g_rho_g(P_0, T_0)   # initial gas density
        eps_x0 = crystal_fraction_eps_x(
            composition, T_0, P_0, InitialConc_H2O, InitialConc_CO2
        )

        eps_m0 = 1 - eps_x0
        rho = rc.rho_m0 * eps_m0 + rc.rho_x0 * eps_x0
        M_tot = V_0 * rho   # Total mass, initial
        M_co2_0 = InitialConc_CO2 * V_0 * rho   # Total mass of CO2, initial
        M_h2o_0 = InitialConc_H2O * V_0 * rho   # Total mass of H2O, initial

        # IC Finder parameters
        param_IC_Finder.Tol = if composition == Silicic()
            1e-9
        else
            1e-8
        end

        eps_g0, X_co20, C_co2, phase = IC_Finder(
            composition, M_h2o_0, M_co2_0, M_tot, P_0, T_0, V_0, rc.rho_m0, param_IC_Finder
        )

        println("IC_Finder done")
        @show eps_g0 X_co20 C_co2
        println("phase: ", phase)
        @info("IC_Finder done: $eps_g0, $X_co20, $C_co2, $phase")
        param_saved_var.phase = phase

        # update initial bulk density (kg/m^3)
        rho_0 = rho_0_f(eps_g0, eps_x0, rho_g0, rc.rho_m0, rc.rho_x0)

        # update solubility
        if phase == 2
            X_co20 = 0.0
        end

        # Calculate the water content (concentration) for inflowing magma contents
        tot_h2o_frac_in = M_h2o_0 / (rho_0 * V_0)   # CHANGE MASS FRACTION FROM XCO2_IN
        tot_co2_frac_in = M_co2_0 / (rho_0 * V_0)
        if param.fluxing
            tot_co2_frac_in =
                param.XCO2_in * param.mm_co2 /
                (param.XCO2_in * param.mm_co2 + (1 - param.XCO2_in) * param.mm_h2o)
            tot_h2o_frac_in = 1 - tot_co2_frac_in
        end
        param.tot_h2o_frac_in = tot_h2o_frac_in
        param.tot_co2_frac_in = tot_co2_frac_in

        tot_Mass_0 = V_0 * rho_0
        tot_Mass_H2O_0 = M_h2o_0
        tot_Mass_CO2_0 = M_co2_0

        if param.single_eruption
            sw.eruption = 1
        end

        # initialize vector to store quantities
        param_saved_var.storeTime = Vector{Float64}([0])
        param_saved_var.storeTemp = Vector{Float64}([T_0])

        @info("sw: $(sw)")
        @info("IC_Finder parameters: $(param_IC_Finder)")
        @info("params: $(param)")

        stopChamber_MT′(out, u, t, int) = stopChamber_MT(out, u, t, int, sw, param)
        affect!′(int, idx) = affect!(int, idx, sw, param, param_saved_var, param_IC_Finder)

        tspan = (0, end_time)
        IC = [
            P_0,
            T_0,
            eps_g0,
            V_0,
            rc.rho_m0,
            rc.rho_x0,
            X_co20,
            tot_Mass_0,
            tot_Mass_H2O_0,
            tot_Mass_CO2_0,
        ]
        println("tspan: ", tspan)
        println("IC: ", IC)
        @info("IC: $IC")
        cb = VectorContinuousCallback(
            stopChamber_MT′, affect!′, 8; rootfind=SciMLBase.RightRootFind
        )
        prob = ODEProblem(odeChamber, IC, tspan, (param, param_saved_var, sw))
        sol = solve(
            prob,
            methods[method];
            callback=cb,
            reltol=odesetting.reltol,
            abstol=odesetting.abstol,
            dt=odesetting.first_step,
            dtmax=odesetting.max_step,
        )
    end
    println(to)
    @info(to)
    reset_timer!(to)
    close(io)
    df = DataFrame(sol)
    write_csv(df, path)
    plot_figs(df, path)

    println(".. Done!")
    return df
end
