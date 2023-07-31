using Chamber
include("./solver_methods.jl")

"""
    chamber(composition::Union{Silicic,Mafic}, end_time::Float64, log_volume_km3::Float64, InitialConc_H2O::Float64, InitialConc_CO2::Float64, log_vfr::Float64, depth::Float64, output_dirname::String; kwargs...)

Simulate the eruption of a volcano using a model for the frequency of eruptions of upper crustal magma chambers based on Degruyter and Huber (2014).

## Arguments
- `composition`: The magma composition. Use `Silicic()` for rhyolite composition (arc magma) or `Mafic()` for basalt composition (rift).
- `end_time`: Maximum magma chamber evolution duration in seconds.
- `log_volume_km3`: The initial volume of the chamber in logarithmic scale. The actual initial chamber volume is calculated as 10^(`log_volume_km3`) in km³.
- `InitialConc_H2O`: The initial weight fraction of water in the magma (exsolved + dissolved).
- `InitialConc_CO2`: The initial weight fraction of CO₂ in the magma (exsolved + dissolved).
- `log_vfr`: Magma recharge rate in km³/yr calculated as 10^(`log_vfr`).
- `depth`: Depth of the magma chamber in meters.
- `output_dirname`(optional): Name of the output directory. Defaults to current timestamp.

## Keyword Arguments
- `plotfig`(optional): (default: `true`). Generate and plot figures for each result if true.

## Returns
A `DataFrame` containing the solution with columns:
- `time`: Simulation timestamps in seconds.
- `P+dP`: Pressure in Pa.
- `T`: Temperature in K.
- `eps_g`: Gas volume fraction.
- `V`: Volume of the magma chamber in m³.
- `rho_m`: Density of the melt in kg/m³.
- `rho_x`: Density of magma crystal in kg/m³.
- `X_CO2`: Mole fraction of CO2 in the gas.
- `total_mass`: Total mass of magma chamber in kg.
- `total_mass_H2O`: Total mass of water in the magma in kg.
- `total_mass_CO2`: Total mass of CO₂ in the magma in kg.

## Outputs
A directory named after `output_dirname` or the default value, containing the following files:
- `out.csv`: a CSV file containing the solution columns listed above.
- `eruptions.csv`, A CSV file containing the datas of eruptions with the following columns: time_of_eruption (sec), duration_of_eruption (sec), mass_erupted (kg) and volume_erupted (km³).
- Figures for P+dP(t), T(t), eps_g(t), V(t), X_CO2(t), total_mass(t).

## References
- W. Degruyter and C. Huber (2014). A model for eruption frequency of upper crustal silicic magma chambers. Earth Planet. Sci. Lett. https://doi.org/10.1016/j.epsl.2014.06.047

## Examples
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

# Run a simulation with mafic magma chamber, with custom directory name "MyDirname"
julia> composition = Mafic()
julia> end_time = 3e9
julia> log_volume_km3 = 0.2
julia> InitialConc_H2O = 0.01
julia> InitialConc_CO2 = 0.001
julia> log_vfr = -3.3
julia> depth = 8e3
julia> output_dirname = "MyDirname"

julia> chamber(composition, end_time, log_volume_km3, InitialConc_H2O, InitialConc_CO2, log_vfr, depth, output_dirname)
```
"""
@memoize function chamber(
    composition::Union{Silicic,Mafic},
    end_time::Float64,
    log_volume_km3::Float64,
    InitialConc_H2O::Float64,
    InitialConc_CO2::Float64,
    log_vfr::Float64,
    depth::Float64,
    output_dirname::String=get_timestamp(),
    ;
    method::String="CVODE_BDF",
    rheol::String="old",
    plotfig::Bool=true,
)::DataFrame
    datetime = get_timestamp()
    composition_str = string(typeof(composition))
    path = joinpath(pwd(), output_dirname)
    mkdir(path)
    println(
        "$(Threads.nthreads() > 1 ? "(thread $(Threads.threadid()) / $(Threads.nthreads())) " : "")Output path: $path",
    )
    io = open("$path/$datetime.log", "w+")
    df = with_logger(SimpleLogger(io)) do
        @info(
            "Arguments:",
            composition,
            end_time,
            log_volume_km3,
            InitialConc_H2O,
            InitialConc_CO2,
            log_vfr,
            depth,
            path
        )

        to = TimerOutput()
        @timeit to "vol$(log_volume_km3)_h2o$(InitialConc_H2O)_gas$(InitialConc_CO2)_vfr$(log_vfr)_dep$(depth)" begin
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

            c = ConstantValues{Float64}()
            param_IC_Finder = ParamICFinder{Float64}()
            param_saved_var = ParamSaved{Float64}()
            sw = SW{Int8}()
            odesetting = OdeSetting{Float64}()
            erupt_saved = EruptSaved{Float64}()

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

            T_0 = find_liq(
                composition, InitialConc_H2O, InitialConc_CO2, P_0, param.ini_eps_x
            )

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
                composition,
                M_h2o_0,
                M_co2_0,
                M_tot,
                P_0,
                T_0,
                V_0,
                rc.rho_m0,
                param_IC_Finder,
            )

            @info("First IC_Finder done: ", eps_g0, X_co20, C_co2, phase)
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

            @info("IC_Finder parameters: $(param_IC_Finder)")

            stopChamber_MT′(out, u, t, int) = stopChamber_MT(out, u, t, int, sw, param)
            function affect!′(int, idx)
                return affect!(
                    int, idx, sw, param, param_saved_var, param_IC_Finder, erupt_saved
                )
            end

            timespan = (0, end_time)
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
            @info("ODE solver settings: ", method, odesetting, IC, timespan, param, sw)
            cb = VectorContinuousCallback(
                stopChamber_MT′, affect!′, 8; rootfind=SciMLBase.RightRootFind
            )
            prob = ODEProblem(odeChamber, IC, timespan, (param, param_saved_var, sw))
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
        print_timer_log(io, to)
        close(io)
        df = DataFrame(sol)
        write_csv(df, erupt_saved, path)
        if plotfig
            plot_figs(df, path)
        end
        return df
    end
    return df
end

"""
    chamber(composition::Union{Silicic,Mafic}, end_time::Float64, log_volume_km3_vector::Union{Float64,Vector{Float64}}, InitialConc_H2O_vector::Union{Float64,Vector{Float64}}, InitialConc_CO2_vector::Union{Float64,Vector{Float64}}, log_vfr_vector::Union{Float64,Vector{Float64}}, depth_vector::Union{Float64,Vector{Float64}}, output_dirname::String; kwargs...)

Simulate the eruption of a volcano using a model for the frequency of eruptions of upper crustal magma chambers based on Degruyter and Huber (2014).

## Arguments
- `composition`: The magma composition. Use `Silicic()` for rhyolite composition (arc magma) or `Mafic()` for basalt composition (rift).
- `end_time`: Maximum magma chamber evolution duration in seconds.
- `log_volume_km3`: The initial volume of the chamber in logarithmic scale. The actual initial chamber volume is calculated as 10^(`log_volume_km3`) in km³.
- `InitialConc_H2O`: The initial weight fraction of water in the magma (exsolved + dissolved).
- `InitialConc_CO2`: The initial weight fraction of CO₂ in the magma (exsolved + dissolved).
- `log_vfr`: Magma recharge rate in km³/yr calculated as 10^(`log_vfr`).
- `depth`: Depth of the magma chamber in meters.
- `output_dirname`(optional): Name of the output directory. Defaults to current timestamp.

## Keyword Arguments
- `plotfig`(optional): (default: `true`). Generate and plot figures for each result if true.

## Returns
A `DataFrame` containing the solution with columns:
- `time`: Simulation timestamps in seconds.
- `P+dP`: Pressure in Pa.
- `T`: Temperature in K.
- `eps_g`: Gas volume fraction.
- `V`: Volume of the magma chamber in m³.
- `rho_m`: Density of the melt in kg/m³.
- `rho_x`: Density of magma crystal in kg/m³.
- `X_CO2`: Mole fraction of CO2 in the gas.
- `total_mass`: Total mass of magma chamber in kg.
- `total_mass_H2O`: Total mass of water in the magma in kg.
- `total_mass_CO2`: Total mass of CO₂ in the magma in kg.

## Outputs
A directory named after `output_dirname` or the default value, containing the following files:
- `out.csv`: a CSV file containing the solution columns listed above.
- `eruptions.csv`, A CSV file containing the datas of eruptions with the following columns: time_of_eruption (sec), duration_of_eruption (sec), mass_erupted (kg) and volume_erupted (km³).
- Figures for P+dP(t), T(t), eps_g(t), V(t), X_CO2(t), total_mass(t).

## References
- W. Degruyter and C. Huber (2014). A model for eruption frequency of upper crustal silicic magma chambers. Earth Planet. Sci. Lett. https://doi.org/10.1016/j.epsl.2014.06.047

## Examples
```
# Run a simulation with silicic magma chamber
julia> composition = Silicic()
julia> end_time = 3e9
julia> log_volume_km3 = 0.2
julia> InitialConc_H2O = [0.04, 0.045]
julia> InitialConc_CO2 = [0.001, 0.0012]
julia> log_vfr = -3.3
julia> depth = 8e3

julia> chamber(composition, end_time, log_volume_km3, InitialConc_H2O, InitialConc_CO2, log_vfr, depth)

# Run a simulation with mafic magma chamber, with custom directory name "MyDirname"
julia> composition = Mafic()
julia> end_time = 3e9
julia> log_volume_km3 = 0.2
julia> InitialConc_H2O = [0.01, 0.012]
julia> InitialConc_CO2 = 0.001
julia> log_vfr = -3.3
julia> depth = [7e3, 8e3]
julia> output_dirname = "MyDirname"

julia> chamber(composition, end_time, log_volume_km3, InitialConc_H2O, InitialConc_CO2, log_vfr, depth, output_dirname)
```
"""
function chamber(
    composition::Union{Silicic,Mafic},
    end_time::Float64,
    log_volume_km3_vector::Union{Float64,Vector{Float64}},
    InitialConc_H2O_vector::Union{Float64,Vector{Float64}},
    InitialConc_CO2_vector::Union{Float64,Vector{Float64}},
    log_vfr_vector::Union{Float64,Vector{Float64}},
    depth_vector::Union{Float64,Vector{Float64}},
    output_dirname::String=get_timestamp(),
    ;
    plotfig::Bool=true,
)::String
    check_for_duplicates(
        log_volume_km3_vector,
        InitialConc_H2O_vector,
        InitialConc_CO2_vector,
        log_vfr_vector,
        depth_vector,
    )
    path0 = joinpath(pwd(), output_dirname)
    mkdir(path0)
    to = TimerOutput()
    df_outputs = Vector{ChamberOutput}(
        undef,
        length(log_volume_km3_vector) *
        length(InitialConc_H2O_vector) *
        length(InitialConc_CO2_vector) *
        length(log_vfr_vector) *
        length(depth_vector),
    )
    @timeit to "chamber" begin
        Threads.@threads for (
            idx, (log_volume_km3, InitialConc_H2O, InitialConc_CO2, log_vfr, depth)
        ) in collect(
            enumerate(
                Iterators.product(
                    log_volume_km3_vector,
                    InitialConc_H2O_vector,
                    InitialConc_CO2_vector,
                    log_vfr_vector,
                    depth_vector,
                ),
            ),
        )
            dataset = "vol$(log_volume_km3)_h2o$(InitialConc_H2O)_gas$(InitialConc_CO2)_vfr$(log_vfr)_dep$(depth)"
            df = chamber(
                composition,
                end_time,
                log_volume_km3,
                InitialConc_H2O,
                InitialConc_CO2,
                log_vfr,
                depth,
                joinpath(output_dirname, dataset);
                plotfig=false,
            )
            df_outputs[idx] = ChamberOutput(df, joinpath(output_dirname, dataset))
        end
        if plotfig
            for output in df_outputs
                plot_figs(output.df, output.path)
            end
        end
    end
    io0 = open("$path0/$output_dirname.log", "w+")
    with_logger(SimpleLogger(io0)) do
        @info(
            "Arguments:",
            composition,
            end_time,
            log_volume_km3_vector,
            InitialConc_H2O_vector,
            InitialConc_CO2_vector,
            log_vfr_vector,
            depth_vector,
            path0
        )
    end
    print_timer_log(io0, to)
    close(io0)
    return output_dirname
end
