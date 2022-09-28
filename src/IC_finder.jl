module IC_finder
using Roots
include("initial.jl")
using .initial: eos_g, crystal_fraction_silicic, exsolve_silicic, exsolve3_silicic, crystal_fraction_mafic, exsolve_mafic, exsolve3_mafic

function mco2_dissolved_sat_silicic(X, P, T)
    P_MPa = P/1e6
    # function & coefficients from Liu et al 2005
    Pc = P_MPa*X
    Pw = P_MPa*(1-X)
    c1 = 5668
    c2 = -55.99
    c3 = 0.4133
    c4 = 2.041e-3
    sol = real(Pc*(c1+c2*Pw)/T+Pc*(c3*Complex(Pw)^0.5+c4*Complex(Pw)^1.5))
    sol = sol/1e6
    return sol
end

function mco2_dissolved_sat_mafic(X, P, T)
    P_MPa = P/1e6
    # function & coefficients from Liu et al 2005
    Pc = P_MPa*X
    Pw = P_MPa*(1-X)
    c1 = 5668
    c2 = -55.99
    c3 = 0.4133
    c4 = 2.041e-3
    sol = real(Pc*(c1+c2*Pw)/T+Pc*(c3*Complex(Pw)^0.5+c4*Complex(Pw)^1.5))
    sol = sol/1e6
    return sol
end

function meq_water_silicic(X, P, T)
    P = P/1e6
    Pw = (1-X)*P
    Pc = X*P
    T_C = T
    a1 = 354.94
    a2 = 9.623
    a3 = -1.5223
    a4 = 1.2439e-3
    a5 = -1.084e-4
    a6 = -1.362e-5
    sol = real((a1*Complex(Pw)^0.5+a2*Pw+a3*Complex(Pw)^1.5)/T_C+a4*Complex(Pw)^1.5+Pc*(a5*Complex(Pw)^0.5+a6*Pw))
    sol = sol/100
    return sol
end

function meq_water_mafic(X, P, T)
    P = P/1e6
    T_C = T-273.15
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
    sol = b1 + b2*T_C + b3*X + b4*P + b5*T_C*X + b6*T_C*P + b7*X*P + b8*T_C^2 + b9*X^2 + b10*P^2
    sol = sol/100
    return sol
end

function IC_Finder_silicic(M_h2o, M_co2, M_tot, P, T, V, rho_m, mm_co2, mm_h2o, param_IC)
    to = get_timer("share")
    ## IC Finder parameters
    max_count = param_IC["max_count"]
    Tol = param_IC["Tol"]
    min_eps_g = param_IC["min_eps_g"]
    eps_g_guess_ini = param_IC["eps_g_guess_ini"]
    X_co2_guess_ini = param_IC["X_co2_guess_ini"]
    fraction = param_IC["fraction"]
    delta_X_co2 = param_IC["delta_X_co2"]
    ## ------------------------------

    rho_g = eos_g(P, T)["rho_g"]
    m_h2o_tot = M_h2o/M_tot
    m_co2_tot = M_co2/M_tot
    @timeit to "crystal_fraction_silicic" eps_x0 = crystal_fraction_silicic(T, P, m_h2o_tot, m_co2_tot)[1]

    # Fixing total mass of volatiles set in MainChamber in real model
    eps_m0 = 1 - eps_x0
    m_h2o_melt = M_h2o/(V*rho_m*eps_m0)
    m_co2_melt = M_co2/(V*rho_m*eps_m0)

    # CHECK IF SATURATED
    @timeit to "exsolve_silicic" m_eq_max = exsolve_silicic(P, T, 0.0)[1]

    if m_h2o_melt > m_eq_max
        phase = 3
    else
        @timeit to "exsolve3_silicic" C_co2_sat = exsolve3_silicic(P, T, m_h2o_melt)[1]
        if m_co2_melt > C_co2_sat  ###
            phase = 3
        else
            phase = 2
        end
    end

    if phase == 2
        X_co20 = 0
        eps_g0 = 0
        mco2_diss = m_co2_melt

    else
        # First Guesses
        eps_g_guess = eps_g_guess_ini
        X_co2_guess = X_co2_guess_ini
        eps_g0 = eps_g_guess
        X_co20 = X_co2_guess
        X_co2_prev = X_co20 + 2*Tol
        eps_g_prev = eps_g0 + 2*Tol
        count = 0

        while ((abs((X_co20-X_co2_prev)/X_co2_prev)>Tol || abs((eps_g0-eps_g_prev)/eps_g_prev)>Tol) && count<max_count)
            @timeit to "crystal_fraction_silicic" eps_x0 = crystal_fraction_silicic(T, P, m_h2o_tot, m_co2_tot)[1]
            fun(x) = real(mco2_dissolved_sat_silicic(x, P, T)*(1-eps_g0-eps_x0)*V*rho_m+eps_g0*rho_g*V*x*mm_co2/(x*mm_co2+(1-x)*mm_h2o)-M_co2)
            X_co2_prev = X_co20
            fx = ZeroProblem(fun, X_co2_prev)
            X_co20 = solve(fx, xatol=Tol, atol=Tol, maxiters=100)

            if X_co20 < 0 || X_co20 > 1
                X_co20 = -1
                break
            end

            Xmean = (1-fraction)*X_co2_prev + fraction*X_co20
            X_co20 = Xmean
            @timeit to "meq_water_silicic" mwater_dissolved = meq_water_silicic(X_co20, P, T)
            @timeit to "mco2_dissolved_sat_silicic" mco2_diss = mco2_dissolved_sat_silicic(X_co20, P, T)
            eps_m = 1-eps_g0-eps_x0
            Num = M_co2 - mco2_diss*eps_m*V*rho_m + M_h2o - mwater_dissolved*eps_m*V*rho_m
            Den = rho_g*V
            eps_g_prev = eps_g0
            eps_g0 = (1-fraction)*eps_g_prev + fraction*Num/Den
            count = count + 1
        end
        while (isnan(X_co20) || isnan(eps_g0) || X_co20<0 || eps_g0<0) && eps_g_guess > min_eps_g
            eps_g_guess = eps_g_guess/1.5
            X_co2_guess = X_co2_guess_ini
            X_co2_prev = X_co2_guess+2*Tol
            eps_g_prev = eps_g_guess+2*Tol
            Err_eps_g = 0
            Err_Xco2 = 0
            while (isnan(X_co20) || isnan(eps_g0) || X_co20<0 || eps_g0<0 || Err_eps_g > Tol || Err_Xco2 > Tol) && X_co2_guess <=1
                eps_g0 = eps_g_guess
                X_co20 = X_co2_guess
                count = 0
                while (abs((X_co20-X_co2_prev)/X_co2_prev)>Tol || abs((eps_g0-eps_g_prev)/eps_g_prev)>Tol) && count<max_count || X_co20 > 1
                    @timeit to "crystal_fraction_silicic" eps_x0 = crystal_fraction_silicic(T, P, m_h2o_tot, m_co2_tot)[1]
                    @timeit to "meq_water_silicic" Mass_saturated = meq_water_silicic(X_co20, P, T)*rho_m*V*(1-eps_g0-eps_x0)
                    fun(x) = real(mco2_dissolved_sat_silicic(x, P, T)*(1-eps_g0-eps_x0)*V*rho_m+eps_g0*rho_g*V*x*mm_co2/(x*mm_co2+(1-x)*mm_h2o)-M_co2)
                    X_co2_prev = X_co20
                    fx = ZeroProblem(fun, X_co2_prev)
                    X_co20 = solve(fx, xatol=Tol, atol=Tol, maxiters=100)
                    Xmean = (1-fraction)*X_co2_prev + fraction*X_co20
                    X_co20 = Xmean
                    @timeit to "meq_water_silicic" mwater_dissolved = meq_water_silicic(X_co20, P, T)
                    @timeit to "mco2_dissolved_sat_silicic" mco2_diss = mco2_dissolved_sat_silicic(X_co20, P, T)
                    eps_m = 1-eps_g0-eps_x0
                    Num = M_co2-mco2_diss*eps_m*V*rho_m+M_h2o-mwater_dissolved*eps_m*V*rho_m
                    Den = rho_g*V
                    eps_g_prev = eps_g0
                    eps_g0 = (1-fraction)*eps_g_prev + fraction*Num/Den
                    count = count+1
                end
                Err_eps_g = abs((eps_g0-eps_g_prev)/eps_g_prev)
                Err_Xco2 = abs((X_co20-X_co2_prev)/X_co2_prev)
                X_co2_guess = X_co2_guess+delta_X_co2
            end
        end
        println("  IC_Finder  eps_g0: $eps_g0, X_co20: $X_co20, count: $count")
        if eps_g0 <= 0 || X_co20 < 0 || real(eps_g0) != eps_g0 || real(X_co20) != X_co20  #if complex number
            X_co20 = 0
            eps_g0 = 0
            mco2_diss = m_co2_melt
            phase = 2
        end
        eps_g0 = real(eps_g0)
    end
    return [eps_g0, X_co20, mco2_diss, phase]
end

function IC_Finder_mafic(M_h2o_0, M_co2_0, M_tot, P_0, T_0, V_0, rho_m0, mm_co2, mm_h2o, param_IC)
    to = get_timer("share")

    ## IC Finder parameters
    max_count = param_IC["max_count"]
    Tol = param_IC["Tol"]
    min_eps_g = param_IC["min_eps_g"]
    eps_g_guess_ini = param_IC["eps_g_guess_ini"]
    X_co2_guess_ini = param_IC["X_co2_guess_ini"]
    fraction = param_IC["fraction"]
    delta_X_co2 = param_IC["delta_X_co2"]
    ## ------------------------------
    count_fzeros = 0
    rho_g0 = eos_g(P_0, T_0)["rho_g"]

    mH2O = M_h2o_0/M_tot
    mCO2 = M_co2_0/M_tot
    @timeit to "crystal_fraction_mafic" eps_x0 = crystal_fraction_mafic(T_0, P_0, mH2O, mCO2)[1]

    eps_m0 = 1 - eps_x0

    Conc_Water = M_h2o_0/(V_0*rho_m0*eps_m0)
    Conc_co2 = M_co2_0/(V_0*rho_m0*eps_m0)

    # CHECK IF SATURATED
    @timeit to "exsolve_mafic" m_eq_max = exsolve_mafic(P_0, T_0, 0)[1]

    if Conc_Water > m_eq_max
        phase = 3
    else
        @timeit to "exsolve3_mafic" C_co2_sat = exsolve3_mafic(P_0, T_0, Conc_Water)[1]
        if Conc_co2 >  C_co2_sat
            phase = 3
        else
            phase = 2
        end
    end

    if phase == 2
        X_co20 = 0
        eps_g0 = 0
        mco2_diss = Conc_co2
    else
        # First Guesses
        eps_g_guess = eps_g_guess_ini
        X_co2_guess = X_co2_guess_ini
        eps_g0 = eps_g_guess
        X_co20 = X_co2_guess
        X_co2_prev = X_co20 + 2*Tol
        eps_g_prev = eps_g0 + 2*Tol
        count = 0

        while ((abs((X_co20-X_co2_prev)/X_co2_prev)>Tol || abs((eps_g0-eps_g_prev)/eps_g_prev)>Tol) && count<max_count)
            @timeit to "crystal_fraction_mafic" eps_x0 = crystal_fraction_mafic(T_0, P_0, mH2O, mCO2)[1]
            fun(x) = real(mco2_dissolved_sat_mafic(x, P_0, T_0)*(1-eps_g0-eps_x0)*V_0*rho_m0+eps_g0*rho_g0*V_0*x*mm_co2/(x*mm_co2+(1-x)*mm_h2o)-M_co2_0)
            if ~isreal(X_co20) || isnan(X_co20)
                X_co2_prev = 0
            else
                X_co2_prev = X_co20
            end
            fx = ZeroProblem(fun, X_co2_prev)
            X_co20 = solve(fx, xatol=Tol, atol=Tol, maxiters=100)
            count_fzeros = count_fzeros+1
            Xmean = (1-fraction)*X_co2_prev + fraction*X_co20
            X_co20 = Xmean
            @timeit to "meq_water_mafic" mwater_dissolved = meq_water_mafic(X_co20, P_0, T_0)
            @timeit to "mco2_dissolved_sat_mafic" mco2_diss = mco2_dissolved_sat_mafic(X_co20, P_0, T_0)
            eps_m = 1-eps_g0-eps_x0
            Num = M_co2_0 - mco2_diss*eps_m*V_0*rho_m0 + M_h2o_0 - mwater_dissolved*eps_m*V_0*rho_m0
            Den = rho_g0*V_0
            eps_g_prev = eps_g0
            eps_g0 = (1-fraction)*eps_g_prev + fraction*Num/Den
            count = count+1
        end

        Err_eps_g = abs((eps_g0-eps_g_prev)/eps_g_prev)
        Err_Xco2 = abs((X_co20-X_co2_prev)/X_co2_prev)
        while (isnan(X_co20) || isnan(eps_g0) || X_co20<0 || eps_g0<0 || X_co20>1 || Err_eps_g > Tol || Err_Xco2 > Tol) && eps_g_guess > min_eps_g
            eps_g_guess = eps_g_guess/2
            X_co2_guess = X_co2_guess_ini
            X_co2_prev = X_co2_guess+2*Tol
            eps_g_prev = eps_g_guess+2*Tol

            while (isnan(X_co20) || isnan(eps_g0) || X_co20<0 || eps_g0<0 || Err_eps_g > Tol || Err_Xco2 > Tol) && X_co2_guess<=1-delta_X_co2
                eps_g0 = eps_g_guess
                X_co20 = X_co2_guess
                count = 0
                count_loop = 0
                while (abs((X_co20-X_co2_prev)/X_co2_prev)>Tol || abs((eps_g0-eps_g_prev)/eps_g_prev)>Tol) && count<max_count && X_co20<=1
                    count_loop = count_loop+1
                    @timeit to "crystal_fraction_mafic" eps_x0 = crystal_fraction_mafic(T_0, P_0, mH2O, mCO2)[1]
                    fun(x) = real(mco2_dissolved_sat_mafic(x, P_0, T_0)*(1-eps_g0-eps_x0)*V_0*rho_m0+eps_g0*rho_g0*V_0*x*mm_co2/(x*mm_co2+(1-x)*mm_h2o)-M_co2_0)

                    if ~isreal(X_co20) || isnan(X_co20) || ~isreal(eps_g0) || isnan(eps_g0)
                        X_co2_prev = 0
                    else
                        X_co2_prev = X_co20
                    end
                    fx = ZeroProblem(fun, X_co2_prev)
                    X_co20 = solve(fx, xatol=Tol, atol=Tol, maxiters=100)
                    count_fzeros = count_fzeros+1
                    Xmean = (1-fraction)*X_co2_prev + fraction*X_co20
                    X_co20 = Xmean
                    @timeit to "meq_water_mafic" mwater_dissolved = meq_water_mafic(X_co20, P_0, T_0)
                    @timeit to "mco2_dissolved_sat_mafic" mco2_diss = mco2_dissolved_sat_mafic(X_co20, P_0, T_0)

                    eps_m = 1-eps_g0-eps_x0
                    Num = M_co2_0 - mco2_diss*eps_m*V_0*rho_m0+M_h2o_0 - mwater_dissolved*eps_m*V_0*rho_m0
                    Den = rho_g0*V_0
                    eps_g_prev = eps_g0
                    eps_g0 = (1-fraction)*eps_g_prev + fraction*Num/Den
                    count = count+1
                end

                if eps_g_prev > 0
                    Err_eps_g = abs((eps_g0-eps_g_prev)/eps_g_prev)
                else
                    Err_eps_g = abs((eps_g0-eps_g_prev))
                end

                if X_co2_prev > 0
                    Err_Xco2 = abs((X_co20-X_co2_prev)/X_co2_prev)
                else
                    Err_Xco2 = abs((X_co20-X_co2_prev))
                end
                X_co2_guess = X_co2_guess + delta_X_co2
            end
        end

        if eps_g0 <= 0 || X_co20 < 0
            X_co20 = 0
            eps_g0 = 0
            mco2_diss = Conc_co2
            phase = 2
        end
        if count == max_count && (Err_eps_g > Tol || Err_Xco2 > Tol)
            X_co20 = 0
            eps_g0 = 0
            mco2_diss = Conc_co2
            phase = 2
        end
        eps_g0 = real(eps_g0)
    end

    if X_co20 > 1
        X_co20 = 0
        eps_g0 = 0
        phase = 2
    end
    return [eps_g0, X_co20, mco2_diss, phase]
end

end