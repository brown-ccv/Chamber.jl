"""
    IC_Finder(composition::Silicic, M_h2o::Float64, M_co2::Float64, M_tot::Float64, P::Float64, T::Float64, V::Float64, rho_m::Float64, param_IC::ParamICFinder{Float64})::NamedTuple{(:eps_g0, :X_co20, :mco2_diss, :phase),NTuple{4,Float64}}

# Arguments
- `M_h2o`: total mass of water in magma
- `M_co2`: total mass of CO2 in magma
- `M_tot`: total mass of magma
- `P`: Pressure (Pa)
- `T`: Temperature (K)
- `V`: chamber volume (m^3)
- `rho_m`: density of the melt
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
    param_IC::ParamICFinder{Float64},
)::NamedTuple{(:eps_g0, :X_co20, :mco2_diss, :phase),NTuple{4,Float64}}
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

    phase, m_co2_melt = get_phase(composition, P, T, V, rho_m, M_h2o, M_co2, eps_x0)
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
        solve_X_co2′(eps_g0, X_co2_prev) = solve_X_co2(eps_g0, X_co2_prev, P, T, eps_x0, V, rho_m, rho_g, M_co2, Tol)
        while ((abs((X_co20 - X_co2_prev) / X_co2_prev) > Tol || abs((eps_g0 - eps_g_prev) / eps_g_prev) > Tol) && count < max_count)
            X_co2_prev = X_co20
            X_co20 = solve_X_co2′(eps_g0, X_co2_prev)

            mwater_dissolved = meq_water(composition, X_co20, P, T)
            mco2_diss = mco2_dissolved_sat(X_co20, P, T)
            eps_m = 1 - eps_g0 - eps_x0
            Num =
                M_co2 - mco2_diss * eps_m * V * rho_m + M_h2o -
                mwater_dissolved * eps_m * V * rho_m
            Den = rho_g * V
            eps_g_prev = eps_g0
            eps_g0 = (1 - fraction) * eps_g_prev + fraction * Num / Den
            count = count + 1
        end
        while (isnan(X_co20) || isnan(eps_g0) || X_co20 < 0 || eps_g0 < 0) && eps_g_guess > min_eps_g
            eps_g_guess = eps_g_guess / 1.5
            X_co2_guess = X_co2_guess_ini
            X_co2_prev = X_co2_guess + 2 * Tol
            eps_g_prev = eps_g_guess + 2 * Tol
            Err_eps_g = 0.0
            Err_Xco2 = 0.0
            while (isnan(X_co20) || isnan(eps_g0) || X_co20 < 0 || eps_g0 < 0 || Err_eps_g > Tol || Err_Xco2 > Tol) && X_co2_guess <= 1
                eps_g0 = eps_g_guess
                X_co20 = X_co2_guess
                count = 0
                while (abs((X_co20 - X_co2_prev) / X_co2_prev) > Tol ||abs((eps_g0 - eps_g_prev) / eps_g_prev) > Tol) && count < max_count || X_co20 > 1
                    X_co2_prev = X_co20
                    X_co20 = solve_X_co2′(eps_g0, X_co2_prev)

                    mwater_dissolved = meq_water(composition, X_co20, P, T)
                    mco2_diss = mco2_dissolved_sat(X_co20, P, T)
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
        if eps_g0 <= 0 || X_co20 < 0 || real(eps_g0) != eps_g0 || real(X_co20) != X_co20  #avoid complex number
            X_co20 = 0.0
            eps_g0 = 0.0
            mco2_diss = m_co2_melt
            phase = 2.0
        end
        eps_g0 = real(eps_g0)
    end
    return (; eps_g0, X_co20, mco2_diss, phase)
end

"""
    IC_Finder(composition::Mafic, M_h2o::Float64, M_co2::Float64, M_tot::Float64, P::Float64, T::Float64, V::Float64, rho_m::Float64, param_IC::ParamICFinder{Float64})::NamedTuple{(:eps_g0, :X_co20, :mco2_diss, :phase),NTuple{4,Float64}}

# Arguments
- `M_h2o`: total mass of water in magma
- `M_co2`: total mass of CO2 in magma
- `M_tot`: total mass of magma
- `P`: Pressure (Pa)
- `T`: Temperature (K)
- `V`: chamber volume (m^3)
- `rho_m`: density of the melt
"""
function IC_Finder(
    composition::Mafic,
    M_h2o::Float64,
    M_co2::Float64,
    M_tot::Float64,
    P::Float64,
    T::Float64,
    V::Float64,
    rho_m::Float64,
    param_IC::ParamICFinder{Float64},
)::NamedTuple{(:eps_g0, :X_co20, :mco2_diss, :phase),NTuple{4,Float64}}
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
    rho_g = eos_g_rho_g(P, T)

    m_h2o_tot = M_h2o / M_tot
    m_co2_tot = M_co2 / M_tot
    eps_x0 = crystal_fraction_eps_x(composition, T, P, m_h2o_tot, m_co2_tot)

    phase, m_co2_melt = get_phase(composition, P, T, V, rho_m, M_h2o, M_co2, eps_x0)

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
        c = ConstantValues()
        mm_co2, mm_h2o = c.mm_co2, c.mm_h2o
        while (
            (
                abs((X_co20 - X_co2_prev) / X_co2_prev) > Tol ||
                abs((eps_g0 - eps_g_prev) / eps_g_prev) > Tol
            ) && count < max_count
        )
            function fun(x)
                return real(
                    mco2_dissolved_sat(x, P, T) * (1 - eps_g0 - eps_x0) * V * rho_m +
                    eps_g0 * rho_g * V * x * mm_co2 / (x * mm_co2 + (1 - x) * mm_h2o) -
                    M_co2,
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
            mwater_dissolved = meq_water(composition, X_co20, P, T)
            mco2_diss = mco2_dissolved_sat(X_co20, P, T)
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
                    function fun(x)
                        return real(
                            mco2_dissolved_sat(x, P, T) *
                            (1 - eps_g0 - eps_x0) *
                            V *
                            rho_m +
                            eps_g0 * rho_g * V * x * mm_co2 /
                            (x * mm_co2 + (1 - x) * mm_h2o) - M_co2,
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
                    mwater_dissolved = meq_water(composition, X_co20, P, T)
                    mco2_diss = mco2_dissolved_sat(X_co20, P, T)

                    eps_m = 1 - eps_g0 - eps_x0
                    Num =
                        M_co2 - mco2_diss * eps_m * V * rho_m + M_h2o -
                        mwater_dissolved * eps_m * V * rho_m
                    Den = rho_g * V
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
    end

    if X_co20 > 1
        X_co20 = 0.0
        eps_g0 = 0.0
        phase = 2
    end
    return (; eps_g0, X_co20, mco2_diss, phase)
end
