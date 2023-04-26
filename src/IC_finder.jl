"""
    IC_Finder(composition::Silicic, M_h2o::Float64, M_co2::Float64, M_tot::Float64, P::Float64, T::Float64, V::Float64, rho_m::Float64, param_IC::ParamICFinder{Float64})::NamedTuple{(:eps_g0, :X_co20, :mco2_diss, :phase),NTuple{4,Float64}}

An iterative function that uses thermodynamic modeling to determine the gas phase composition of a silicic magma chamber.

## Arguments
- `composition`: `Silicic()`
- `M_h2o`: total mass of water in magma (kg)
- `M_co2`: total mass of CO2 in magma (kg)
- `M_tot`: total mass of magma (kg)
- `P`: Pressure (Pa)
- `T`: Temperature (K)
- `V`: chamber volume (m³)
- `rho_m`: density of the melt (kg/m³)
- `param_IC`: An instance of the `ParamICFinder` composite type containing the parameters for IC_Finder

## Returns
A NamedTuple with the following fields:
- `eps_g0`: Initial gas fraction of the magma chamber
- `X_co20`: Initial mole fraction of CO2 in the gas phase
- `mco2_diss`: Total mass of CO2 dissolved in the melt (kg)
- `phase`: Integer representing the state of the magma

## Details
The function uses a thermodynamic model to find the initial gas fraction (`eps_g0`) and initial mole fraction of CO2 in the gas phase (`X_co20`) for a silicic magma chamber. The model assumes that the chamber contains two phases: a liquid melt and a gas phase. The model calculates the saturation state of the magma with respect to CO2, and determines whether the gas phase contains only CO2 or a mixture of CO2 and H2O. If the chamber is in two-phase state, the model returns `eps_g0 = 0.0`, `X_co20 = 0.0`, `mco2_diss = Total mass of CO2 dissolved in the melt`, and `phase = 2`. If the chamber is in gas or liquid state, the function iteratively solves for `eps_g0` and `X_co20` using the `solve_X_co2` and `get_eps_g` helper functions until the relative difference between `eps_g0` and `X_co20` in two consecutive iterations is less than the tolerance `Tol`, or until the maximum number of iterations `max_count` is reached. 
"""
function IC_Finder(
    composition::Silicic,
    M_h2o::T,
    M_co2::T,
    M_tot::T,
    P::T,
    Temp::T,
    V::T,
    rho_m::T,
    param_IC::ParamICFinder{T},
)::NamedTuple{(:eps_g0, :X_co20, :mco2_diss, :phase),NTuple{4,T}} where {T<:Float64}
    ## IC Finder parameters
    max_count = param_IC.max_count
    Tol = param_IC.Tol
    min_eps_g = param_IC.min_eps_g
    eps_g_guess_ini = param_IC.eps_g_guess_ini
    X_co2_guess_ini = param_IC.X_co2_guess_ini
    delta_X_co2 = param_IC.delta_X_co2

    rho_g = eos_g_rho_g(P, Temp)
    m_h2o_tot, m_co2_tot = M_h2o / M_tot, M_co2 / M_tot
    eps_x0 = crystal_fraction_eps_x(composition, Temp, P, m_h2o_tot, m_co2_tot)
    phase, m_co2_melt = get_phase(composition, P, Temp, V, rho_m, M_h2o, M_co2, eps_x0)

    if phase == 2
        return (eps_g0=0.0, X_co20=0.0, mco2_diss=m_co2_melt, phase=phase)
    end

    # First Guesses
    eps_g_guess = eps_g_guess_ini
    X_co2_guess = X_co2_guess_ini
    eps_g0 = eps_g_guess
    X_co20 = X_co2_guess
    X_co2_prev = X_co20 + 2 * Tol
    eps_g_prev = eps_g0 + 2 * Tol
    count = 0
    function solve_X_co2′(eps_g0, X_co2_prev)
        return solve_X_co2(
            composition, eps_g0, X_co2_prev, P, Temp, eps_x0, V, rho_m, rho_g, M_co2, Tol
        )
    end
    function get_eps_g′(eps_g_prev, X_co20)
        return get_eps_g(
            composition, eps_g_prev, X_co20, P, Temp, eps_x0, V, rho_m, rho_g, M_h2o, M_co2
        )
    end
    while (
        abs((X_co20 - X_co2_prev) / X_co2_prev) > Tol ||
        abs((eps_g0 - eps_g_prev) / eps_g_prev) > Tol
    ) && count < max_count
        X_co2_prev = X_co20
        X_co20 = solve_X_co2′(eps_g0, X_co2_prev)
        if X_co20 == -1
            break
        end
        eps_g_prev = eps_g0
        eps_g0, mco2_diss = get_eps_g′(eps_g_prev, X_co20)
        count += 1
    end
    while (isnan(X_co20) || isnan(eps_g0) || X_co20 < 0 || eps_g0 < 0) &&
        eps_g_guess > min_eps_g
        eps_g_guess = eps_g_guess / 1.5
        X_co2_guess = X_co2_guess_ini
        X_co2_prev = X_co2_guess + 2 * Tol
        eps_g_prev = eps_g_guess + 2 * Tol
        Err_eps_g = 0.0
        Err_Xco2 = 0.0
        while (
            isnan(X_co20) ||
            isnan(eps_g0) ||
            X_co20 < 0 ||
            eps_g0 < 0 ||
            Err_eps_g > Tol ||
            Err_Xco2 > Tol
        ) && X_co2_guess <= 1
            eps_g0 = eps_g_guess
            X_co20 = X_co2_guess
            count = 0
            while (
                abs((X_co20 - X_co2_prev) / X_co2_prev) > Tol ||
                abs((eps_g0 - eps_g_prev) / eps_g_prev) > Tol
            ) && count < max_count || X_co20 > 1
                X_co2_prev = X_co20
                X_co20 = solve_X_co2′(eps_g0, X_co2_prev)
                eps_g_prev = eps_g0
                eps_g0, mco2_diss = get_eps_g′(eps_g_prev, X_co20)
                count += 1
            end
            Err_eps_g = abs((eps_g0 - eps_g_prev) / eps_g_prev)
            Err_Xco2 = abs((X_co20 - X_co2_prev) / X_co2_prev)
            X_co2_guess = X_co2_guess + delta_X_co2
        end
    end
    if eps_g0 <= 0 || X_co20 < 0 || !isreal(eps_g0) || !isreal(X_co20)  #avoid complex number
        X_co20 = 0.0
        eps_g0 = 0.0
        mco2_diss = m_co2_melt
        phase = 2.0
    end
    return (; eps_g0, X_co20, mco2_diss, phase)
end

"""
    IC_Finder(composition::Mafic, M_h2o::Float64, M_co2::Float64, M_tot::Float64, P::Float64, T::Float64, V::Float64, rho_m::Float64, param_IC::ParamICFinder{Float64})::NamedTuple{(:eps_g0, :X_co20, :mco2_diss, :phase),NTuple{4,Float64}}

An iterative function that uses thermodynamic modeling to determine the gas phase composition of a mafic magma chamber.

## Arguments
- `composition`: `Mafic()`
- `M_h2o`: total mass of water in magma (kg)
- `M_co2`: total mass of CO2 in magma (kg)
- `M_tot`: total mass of magma (kg)
- `P`: Pressure (Pa)
- `T`: Temperature (K)
- `V`: chamber volume (m³)
- `rho_m`: density of the melt (kg/m³)
- `param_IC`: An instance of the `ParamICFinder` composite type containing the parameters for IC_Finder

## Returns
A NamedTuple with the following fields:
- `eps_g0`: Initial gas fraction of the magma chamber
- `X_co20`: Initial mole fraction of CO2 in the gas phase
- `mco2_diss`: Total mass of CO2 dissolved in the melt (kg)
- `phase`: Integer representing the state of the magma

## Details
The function uses a thermodynamic model to find the initial gas fraction (`eps_g0`) and initial mole fraction of CO2 in the gas phase (`X_co20`) for a mafic magma chamber. The model assumes that the chamber contains two phases: a liquid melt and a gas phase. The model calculates the saturation state of the magma with respect to CO2, and determines whether the gas phase contains only CO2 or a mixture of CO2 and H2O. If the chamber is in two-phase state, the model returns `eps_g0 = 0.0`, `X_co20 = 0.0`, `mco2_diss = Total mass of CO2 dissolved in the melt`, and `phase = 2`. If the chamber is in gas or liquid state, the function iteratively solves for `eps_g0` and `X_co20` using the `solve_X_co2` and `get_eps_g` helper functions until the relative difference between `eps_g0` and `X_co20` in two consecutive iterations is less than the tolerance `Tol`, or until the maximum number of iterations `max_count` is reached. 
"""
function IC_Finder(
    composition::Mafic,
    M_h2o::T,
    M_co2::T,
    M_tot::T,
    P::T,
    Temp::T,
    V::T,
    rho_m::T,
    param_IC::ParamICFinder{T},
)::NamedTuple{(:eps_g0, :X_co20, :mco2_diss, :phase),NTuple{4,T}} where {T<:Float64}
    ## IC Finder parameters
    max_count = param_IC.max_count
    Tol = param_IC.Tol
    min_eps_g = param_IC.min_eps_g
    eps_g_guess_ini = param_IC.eps_g_guess_ini
    X_co2_guess_ini = param_IC.X_co2_guess_ini
    delta_X_co2 = param_IC.delta_X_co2

    rho_g = eos_g_rho_g(P, Temp)
    m_h2o_tot, m_co2_tot = M_h2o / M_tot, M_co2 / M_tot
    eps_x0 = crystal_fraction_eps_x(composition, Temp, P, m_h2o_tot, m_co2_tot)
    phase, m_co2_melt = get_phase(composition, P, Temp, V, rho_m, M_h2o, M_co2, eps_x0)

    if phase == 2
        return (eps_g0=0.0, X_co20=0.0, mco2_diss=m_co2_melt, phase=phase)
    end

    # First Guesses
    eps_g_guess = eps_g_guess_ini
    X_co2_guess = X_co2_guess_ini
    eps_g0 = eps_g_guess
    X_co20 = X_co2_guess
    X_co2_prev = X_co20 + 2 * Tol
    eps_g_prev = eps_g0 + 2 * Tol
    count = 0
    function solve_X_co2′(eps_g0, X_co2_prev)
        return solve_X_co2(
            composition, eps_g0, X_co2_prev, P, Temp, eps_x0, V, rho_m, rho_g, M_co2, Tol
        )
    end
    function get_eps_g′(eps_g_prev, X_co20)
        return get_eps_g(
            composition, eps_g_prev, X_co20, P, Temp, eps_x0, V, rho_m, rho_g, M_h2o, M_co2
        )
    end
    while (
        abs((X_co20 - X_co2_prev) / X_co2_prev) > Tol ||
        abs((eps_g0 - eps_g_prev) / eps_g_prev) > Tol
    ) && count < max_count
        X_co2_prev = isreal(X_co20) && !isnan(X_co20) ? X_co20 : 0.0
        X_co20 = solve_X_co2′(eps_g0, X_co2_prev)
        eps_g_prev = eps_g0
        eps_g0, mco2_diss = get_eps_g′(eps_g_prev, X_co20)
        count += 1
    end

    Err_eps_g = abs((eps_g0 - eps_g_prev) / eps_g_prev)
    Err_Xco2 = abs((X_co20 - X_co2_prev) / X_co2_prev)
    while (
        isnan(X_co20) ||
        isnan(eps_g0) ||
        X_co20 < 0 ||
        eps_g0 < 0 ||
        X_co20 > 1 ||
        Err_eps_g > Tol ||
        Err_Xco2 > Tol
    ) && eps_g_guess > min_eps_g
        eps_g_guess = eps_g_guess / 2
        X_co2_guess = X_co2_guess_ini
        X_co2_prev = X_co2_guess + 2 * Tol
        eps_g_prev = eps_g_guess + 2 * Tol
        while (
            isnan(X_co20) ||
            isnan(eps_g0) ||
            X_co20 < 0 ||
            eps_g0 < 0 ||
            Err_eps_g > Tol ||
            Err_Xco2 > Tol
        ) && X_co2_guess <= 1 - delta_X_co2
            eps_g0 = eps_g_guess
            X_co20 = X_co2_guess
            count = 0
            while (
                          abs((X_co20 - X_co2_prev) / X_co2_prev) > Tol ||
                          abs((eps_g0 - eps_g_prev) / eps_g_prev) > Tol
                      ) &&
                      count < max_count &&
                      X_co20 <= 1
                X_co2_prev =
                    if (any(isreal, [X_co20, eps_g0]) && any(!isnan, [X_co20, eps_g0]))
                        X_co20
                    else
                        0.0
                    end
                X_co20 = solve_X_co2′(eps_g0, X_co2_prev)
                eps_g_prev = eps_g0
                eps_g0, mco2_diss = get_eps_g′(eps_g_prev, X_co20)
                count += 1
            end

            if eps_g_prev > 0
                Err_eps_g = abs((eps_g0 - eps_g_prev) / eps_g_prev)
            else
                Err_eps_g = abs((eps_g0 - eps_g_prev))
            end

            if X_co2_prev > 0
                Err_Xco2 = abs((X_co20 - X_co2_prev) / X_co2_prev)
            else
                Err_Xco2 = abs((X_co20 - X_co2_prev))
            end
            X_co2_guess = X_co2_guess + delta_X_co2
        end
    end

    if eps_g0 <= 0 || X_co20 < 0
        X_co20 = 0.0
        eps_g0 = 0.0
        mco2_diss = m_co2_melt
        phase = 2
    end
    if count == max_count && (Err_eps_g > Tol || Err_Xco2 > Tol)
        X_co20 = 0.0
        eps_g0 = 0.0
        mco2_diss = m_co2_melt
        phase = 2
    end
    eps_g0 = real(eps_g0)

    if X_co20 > 1
        X_co20 = 0.0
        eps_g0 = 0.0
        phase = 2
    end
    return (; eps_g0, X_co20, mco2_diss, phase)
end
