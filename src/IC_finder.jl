# Functions for initial condition
"""
    mco2_dissolved_sat_silicic(X::Float64, P::Float64, T::Float64)::Float64

# Arguments
`X`: mole fraction of CO2 in gas
`P`: Pressure (Pa)
`T`: Temperature (K)
"""
function mco2_dissolved_sat_silicic(X::Float64, P::Float64, T::Float64)::Float64
    P_MPa = P / 1e6
    # function & coefficients from Liu et al 2005
    Pc = P_MPa * X
    Pw = P_MPa * (1 - X)
    @unpack c1, c2, c3, c4 = Co2PartitionCoeff()
    sol = real(Pc * (c1 + c2 * Pw) / T + Pc * (c3 * Complex(Pw)^0.5 + c4 * Complex(Pw)^1.5))
    sol = sol / 1e6
    return sol
end

"""
    mco2_dissolved_sat_mafic(X::Float64, P::Float64, T::Float64)::Float64

# Arguments
`X`: mole fraction of CO2 in gas
`P`: Pressure (Pa)
`T`: Temperature (K)
"""
function mco2_dissolved_sat_mafic(X::Float64, P::Float64, T::Float64)::Float64
    P_MPa = P / 1e6
    # function & coefficients from Liu et al 2005
    Pc = P_MPa * X
    Pw = P_MPa * (1 - X)
    @unpack c1, c2, c3, c4 = Co2PartitionCoeff()
    sol = real(Pc * (c1 + c2 * Pw) / T + Pc * (c3 * Complex(Pw)^0.5 + c4 * Complex(Pw)^1.5))
    sol = sol / 1e6
    return sol
end

"""
    meq_water_silicic(X::Float64, P::Float64, T::Float64)::Float64

# Arguments
`X`: mole fration of H2O in gas
`P`: Pressure (Pa)
`T`: Temperature (K)
"""
function meq_water_silicic(X::Float64, P::Float64, T::Float64)::Float64
    P = P / 1e6
    Pw = (1 - X) * P
    Pc = X * P
    T_C = T
    a1 = 354.94
    a2 = 9.623
    a3 = -1.5223
    a4 = 1.2439e-3
    a5 = -1.084e-4
    a6 = -1.362e-5
    sol = real(
        (a1 * Complex(Pw)^0.5 + a2 * Pw + a3 * Complex(Pw)^1.5) / T_C +
        a4 * Complex(Pw)^1.5 +
        Pc * (a5 * Complex(Pw)^0.5 + a6 * Pw),
    )
    sol = sol / 100
    return sol
end

"""
    meq_water_mafic(X::Float64, P::Float64, T::Float64)::Float64

# Arguments
`X`: mole fration of H2O in gas
`P`: Pressure (Pa)
`T`: Temperature (K)
"""
function meq_water_mafic(X::Float64, P::Float64, T::Float64)::Float64
    P = P / 1e6
    T_C = T - 273.15
    b1 = 2.99622526644026
    b2 = 0.00322422830627781
    b3 = -9.1389095360385
    b4 = 0.0336065247530767
    b5 = 0.00747236662935722
    b6 = -0.0000150329805347769
    b7 = -0.01233608521548
    b8 = -4.14842647942619e-6
    b9 = -0.655454303068124
    b10 = -7.35270395041104e-6
    sol =
        b1 +
        b2 * T_C +
        b3 * X +
        b4 * P +
        b5 * T_C * X +
        b6 * T_C * P +
        b7 * X * P +
        b8 * T_C^2 +
        b9 * X^2 +
        b10 * P^2
    sol = sol / 100
    return sol
end

"""
    IC_Finder(composition::Silicic, M_h2o::Float64, M_co2::Float64, M_tot::Float64, P::Float64, T::Float64, V::Float64, rho_m::Float64, mm_co2::Float64, mm_h2o::Float64, param_IC::ParamICFinder{Float64})::Vector{Float64}

# Arguments
`M_h2o`: total mass of water in magma
`M_co2`: total mass of CO2 in magma
`M_tot`: total mass of magma
`P`: Pressure (Pa)
`T`: Temperature (K)
`V`: chamber volume (m^3)
`rho_m`: density of the melt
`mm_co2`: molecular mass of CO2
`mm_h2o`: molecular mass of H2O
"""
function IC_Finder(
    composition::Silicic,
    M_h2o::Float64,
    M_co2::Float64,
    M_tot::Float64,
    P::Float64,
    T::Float64,
    V::Float64,
    rho_m::Float64,
    mm_co2::Float64,
    mm_h2o::Float64,
    param_IC::ParamICFinder{Float64},
)::Vector{Float64}
    ## IC Finder parameters
    max_count = param_IC.max_count
    Tol = param_IC.Tol
    min_eps_g = param_IC.min_eps_g
    eps_g_guess_ini = param_IC.eps_g_guess_ini
    X_co2_guess_ini = param_IC.X_co2_guess_ini
    fraction = param_IC.fraction
    delta_X_co2 = param_IC.delta_X_co2
    ## ------------------------------

    rho_g = eos_g_rho_g(P, T)
    m_h2o_tot = M_h2o / M_tot
    m_co2_tot = M_co2 / M_tot
    eps_x0 = crystal_fraction_eps_x(composition, T, P, m_h2o_tot, m_co2_tot)

    # Fixing total mass of volatiles set in MainChamber in real model
    eps_m0 = 1 - eps_x0
    m_h2o_melt = M_h2o / (V * rho_m * eps_m0)
    m_co2_melt = M_co2 / (V * rho_m * eps_m0)

    # CHECK IF SATURATED
    m_eq_max = exsolve_meq(composition, P, T, 0.0)

    if m_h2o_melt > m_eq_max
        phase = 3
    else
        C_co2_sat = exsolve3(composition, P, T, m_h2o_melt)[1]
        if m_co2_melt > C_co2_sat  ###
            phase = 3
        else
            phase = 2
        end
    end

    if phase == 2
        X_co20 = 0.0
        eps_g0 = 0.0
        mco2_diss = m_co2_melt

    else
        # First Guesses
        eps_g_guess = eps_g_guess_ini
        X_co2_guess = X_co2_guess_ini
        eps_g0 = eps_g_guess
        X_co20 = X_co2_guess
        X_co2_prev = X_co20 + 2 * Tol
        eps_g_prev = eps_g0 + 2 * Tol
        count = 0

        while (
            (
                abs((X_co20 - X_co2_prev) / X_co2_prev) > Tol ||
                abs((eps_g0 - eps_g_prev) / eps_g_prev) > Tol
            ) && count < max_count
        )
            eps_x0 = crystal_fraction_eps_x(composition, T, P, m_h2o_tot, m_co2_tot)
            function fun(x)
                return real(
                    mco2_dissolved_sat_silicic(x, P, T) *
                    (1 - eps_g0 - eps_x0) *
                    V *
                    rho_m +
                    eps_g0 * rho_g * V * x * mm_co2 / (x * mm_co2 + (1 - x) * mm_h2o) -
                    M_co2,
                )
            end
            X_co2_prev = X_co20
            fx = ZeroProblem(fun, X_co2_prev)
            X_co20 = solve(fx; xatol=Tol, atol=Tol, maxiters=100)

            if X_co20 < 0 || X_co20 > 1
                X_co20 = -1.0
                break
            end

            Xmean = (1 - fraction) * X_co2_prev + fraction * X_co20
            X_co20 = Xmean
            mwater_dissolved = meq_water_silicic(X_co20, P, T)
            mco2_diss = mco2_dissolved_sat_silicic(X_co20, P, T)
            eps_m = 1 - eps_g0 - eps_x0
            Num =
                M_co2 - mco2_diss * eps_m * V * rho_m + M_h2o -
                mwater_dissolved * eps_m * V * rho_m
            Den = rho_g * V
            eps_g_prev = eps_g0
            eps_g0 = (1 - fraction) * eps_g_prev + fraction * Num / Den
            count = count + 1
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
                    eps_x0 = crystal_fraction_eps_x(composition, T, P, m_h2o_tot, m_co2_tot)
                    function fun(x)
                        return real(
                            mco2_dissolved_sat_silicic(x, P, T) *
                            (1 - eps_g0 - eps_x0) *
                            V *
                            rho_m +
                            eps_g0 * rho_g * V * x * mm_co2 /
                            (x * mm_co2 + (1 - x) * mm_h2o) - M_co2,
                        )
                    end
                    X_co2_prev = X_co20
                    fx = ZeroProblem(fun, X_co2_prev)
                    X_co20 = solve(fx; xatol=Tol, atol=Tol, maxiters=100)
                    Xmean = (1 - fraction) * X_co2_prev + fraction * X_co20
                    X_co20 = Xmean
                    mwater_dissolved = meq_water_silicic(X_co20, P, T)
                    mco2_diss = mco2_dissolved_sat_silicic(X_co20, P, T)
                    eps_m = 1 - eps_g0 - eps_x0
                    Num =
                        M_co2 - mco2_diss * eps_m * V * rho_m + M_h2o -
                        mwater_dissolved * eps_m * V * rho_m
                    Den = rho_g * V
                    eps_g_prev = eps_g0
                    eps_g0 = (1 - fraction) * eps_g_prev + fraction * Num / Den
                    count = count + 1
                end
                Err_eps_g = abs((eps_g0 - eps_g_prev) / eps_g_prev)
                Err_Xco2 = abs((X_co20 - X_co2_prev) / X_co2_prev)
                X_co2_guess = X_co2_guess + delta_X_co2
            end
        end
        println("  IC_Finder  eps_g0: $eps_g0, X_co20: $X_co20, count: $count")
        if eps_g0 <= 0 || X_co20 < 0 || real(eps_g0) != eps_g0 || real(X_co20) != X_co20  #if complex number
            X_co20 = 0.0
            eps_g0 = 0.0
            mco2_diss = m_co2_melt
            phase = 2
        end
        eps_g0 = real(eps_g0)
    end
    return [eps_g0, X_co20, mco2_diss, phase]
end

"""
    IC_Finder(composition::Mafic, M_h2o_0::Float64, M_co2_0::Float64, M_tot::Float64, P_0::Float64, T_0::Float64, V_0::Float64, rho_m0::Float64, mm_co2::Float64, mm_h2o::Float64, param_IC::ParamICFinder{Float64})::Vector{Float64}

# Arguments
`M_h2o_0`: total mass of water in magma
`M_co2_0`: total mass of CO2 in magma
`M_tot`: total mass of magma
`P_0`: Pressure (Pa)
`T_0`: Temperature (K)
`V_0`: chamber volume (m^3)
`rho_m0`: density of the melt
`mm_co2`: molecular mass of CO2
`mm_h2o`: molecular mass of H2O
"""
function IC_Finder(
    composition::Mafic,
    M_h2o_0::Float64,
    M_co2_0::Float64,
    M_tot::Float64,
    P_0::Float64,
    T_0::Float64,
    V_0::Float64,
    rho_m0::Float64,
    mm_co2::Float64,
    mm_h2o::Float64,
    param_IC::ParamICFinder{Float64},
)::Vector{Float64}
    ## IC Finder parameters
    max_count = param_IC.max_count
    Tol = param_IC.Tol
    min_eps_g = param_IC.min_eps_g
    eps_g_guess_ini = param_IC.eps_g_guess_ini
    X_co2_guess_ini = param_IC.X_co2_guess_ini
    fraction = param_IC.fraction
    delta_X_co2 = param_IC.delta_X_co2
    ## ------------------------------
    count_fzeros = 0
    rho_g0 = eos_g_rho_g(P_0, T_0)

    mH2O = M_h2o_0 / M_tot
    mCO2 = M_co2_0 / M_tot
    eps_x0 = crystal_fraction_eps_x(composition, T_0, P_0, mH2O, mCO2)

    eps_m0 = 1 - eps_x0

    Conc_Water = M_h2o_0 / (V_0 * rho_m0 * eps_m0)
    Conc_co2 = M_co2_0 / (V_0 * rho_m0 * eps_m0)

    # CHECK IF SATURATED
    m_eq_max = exsolve_meq(composition, P_0, T_0, 0.0)

    if Conc_Water > m_eq_max
        phase = 3
    else
        C_co2_sat = exsolve3(composition, P_0, T_0, Conc_Water)[1]
        if Conc_co2 > C_co2_sat
            phase = 3
        else
            phase = 2
        end
    end

    if phase == 2
        X_co20 = 0.0
        eps_g0 = 0.0
        mco2_diss = Conc_co2
    else
        # First Guesses
        eps_g_guess = eps_g_guess_ini
        X_co2_guess = X_co2_guess_ini
        eps_g0 = eps_g_guess
        X_co20 = X_co2_guess
        X_co2_prev = X_co20 + 2 * Tol
        eps_g_prev = eps_g0 + 2 * Tol
        count = 0

        while (
            (
                abs((X_co20 - X_co2_prev) / X_co2_prev) > Tol ||
                abs((eps_g0 - eps_g_prev) / eps_g_prev) > Tol
            ) && count < max_count
        )
            eps_x0 = crystal_fraction_eps_x(composition, T_0, P_0, mH2O, mCO2)
            function fun(x)
                return real(
                    mco2_dissolved_sat_mafic(x, P_0, T_0) *
                    (1 - eps_g0 - eps_x0) *
                    V_0 *
                    rho_m0 +
                    eps_g0 * rho_g0 * V_0 * x * mm_co2 / (x * mm_co2 + (1 - x) * mm_h2o) -
                    M_co2_0,
                )
            end
            if ~isreal(X_co20) || isnan(X_co20)
                X_co2_prev = 0.0
            else
                X_co2_prev = X_co20
            end
            fx = ZeroProblem(fun, X_co2_prev)
            X_co20 = solve(fx; xatol=Tol, atol=Tol, maxiters=100)
            count_fzeros = count_fzeros + 1
            Xmean = (1 - fraction) * X_co2_prev + fraction * X_co20
            X_co20 = Xmean
            mwater_dissolved = meq_water_mafic(X_co20, P_0, T_0)
            mco2_diss = mco2_dissolved_sat_mafic(X_co20, P_0, T_0)
            eps_m = 1 - eps_g0 - eps_x0
            Num =
                M_co2_0 - mco2_diss * eps_m * V_0 * rho_m0 + M_h2o_0 -
                mwater_dissolved * eps_m * V_0 * rho_m0
            Den = rho_g0 * V_0
            eps_g_prev = eps_g0
            eps_g0 = (1 - fraction) * eps_g_prev + fraction * Num / Den
            count = count + 1
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
                count_loop = 0
                while (
                              abs((X_co20 - X_co2_prev) / X_co2_prev) > Tol ||
                              abs((eps_g0 - eps_g_prev) / eps_g_prev) > Tol
                          ) &&
                          count < max_count &&
                          X_co20 <= 1
                    count_loop = count_loop + 1
                    eps_x0 = crystal_fraction_eps_x(composition, T_0, P_0, mH2O, mCO2)
                    function fun(x)
                        return real(
                            mco2_dissolved_sat_mafic(x, P_0, T_0) *
                            (1 - eps_g0 - eps_x0) *
                            V_0 *
                            rho_m0 +
                            eps_g0 * rho_g0 * V_0 * x * mm_co2 /
                            (x * mm_co2 + (1 - x) * mm_h2o) - M_co2_0,
                        )
                    end

                    if ~isreal(X_co20) || isnan(X_co20) || ~isreal(eps_g0) || isnan(eps_g0)
                        X_co2_prev = 0.0
                    else
                        X_co2_prev = X_co20
                    end
                    fx = ZeroProblem(fun, X_co2_prev)
                    X_co20 = solve(fx; xatol=Tol, atol=Tol, maxiters=100)
                    count_fzeros = count_fzeros + 1
                    Xmean = (1 - fraction) * X_co2_prev + fraction * X_co20
                    X_co20 = Xmean
                    mwater_dissolved = meq_water_mafic(X_co20, P_0, T_0)
                    mco2_diss = mco2_dissolved_sat_mafic(X_co20, P_0, T_0)

                    eps_m = 1 - eps_g0 - eps_x0
                    Num =
                        M_co2_0 - mco2_diss * eps_m * V_0 * rho_m0 + M_h2o_0 -
                        mwater_dissolved * eps_m * V_0 * rho_m0
                    Den = rho_g0 * V_0
                    eps_g_prev = eps_g0
                    eps_g0 = (1 - fraction) * eps_g_prev + fraction * Num / Den
                    count = count + 1
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
            mco2_diss = Conc_co2
            phase = 2
        end
        if count == max_count && (Err_eps_g > Tol || Err_Xco2 > Tol)
            X_co20 = 0.0
            eps_g0 = 0.0
            mco2_diss = Conc_co2
            phase = 2
        end
        eps_g0 = real(eps_g0)
    end

    if X_co20 > 1
        X_co20 = 0.0
        eps_g0 = 0.0
        phase = 2
    end
    return [eps_g0, X_co20, mco2_diss, phase]
end
