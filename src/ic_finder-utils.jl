"""
    mco2_dissolved_sat(X::Float64, P::Float64, T::Float64)::Float64

# Arguments
`X`: mole fraction of CO2 in gas
`P`: Pressure (Pa)
`T`: Temperature (K)
"""
function mco2_dissolved_sat(X::Float64, P::Float64, T::Float64)::Float64
    P_MPa = P / 1e6
    # function & coefficients from Liu et al 2005
    Pc = P_MPa * X
    Pw = P_MPa * (1 - X)
    @unpack c1, c2, c3, c4 = Co2PartitionCoeff()
    sol = real(Pc * (c1 + c2 * Pw) / T + Pc * (c3 * Complex(Pw)^0.5 + c4 * Complex(Pw)^1.5)) / 1e6
    return sol
end

"""
    meq_water(composition::Silicic, X::Float64, P::Float64, T::Float64)::Float64

# Arguments
`X`: mole fration of H2O in gas
`P`: Pressure (Pa)
`T`: Temperature (K)
"""
function meq_water(composition::Silicic, X::Float64, P::Float64, T::Float64)::Float64
    P = P / 1e6
    Pw = (1 - X) * P
    Pc = X * P
    @unpack h1, h2, h3, h4, h5, h6 = Exsolve3Silicic()
    sol = real(
        (h1 * Complex(Pw)^0.5 + h2 * Pw + h3 * Complex(Pw)^1.5) / T +
        h4 * Complex(Pw)^1.5 +
        Pc * (h5 * Complex(Pw)^0.5 + h6 * Pw),
    ) / 100
    return sol
end

"""
    meq_water(composition::Mafic, X::Float64, P::Float64, T::Float64)::Float64

# Arguments
`X`: mole fration of H2O in gas
`P`: Pressure (Pa)
`T`: Temperature (K)
"""
function meq_water(composition::Mafic, X::Float64, P::Float64, T::Float64)::Float64
    P = P / 1e6
    T_C = T - 273.15
    @unpack h1, h2, h3, h4, h5, h6, h7, h8, h9, h10 = Exsolve3Mafic()
    sol =
        h1 +
        h2 * T_C +
        h3 * X +
        h4 * P +
        h5 * T_C * X +
        h6 * T_C * P +
        h7 * X * P +
        h8 * T_C^2 +
        h9 * X^2 +
        h10 * P^2
    return sol / 100
end

"""
    get_phase(composition::Union{Silicic,Mafic}, P::Float64, T::Float64, V::Float64, rho_m::Float64, M_h2o::Float64, M_co2::Float64, eps_x0::Float64)::NamedTuple{(:phase, :m_co2_melt),NTuple{2,Float64}}

- This function is used within the `IC_Finder` function to determine the initial phase.

# Returns
- `phase`: 2 or 3
- `m_co2_melt`
"""
function get_phase(composition::Union{Silicic,Mafic}, P::Float64, T::Float64, V::Float64, rho_m::Float64, M_h2o::Float64, M_co2::Float64, eps_x0::Float64)::NamedTuple{(:phase, :m_co2_melt),NTuple{2,Float64}}
    # Fixing total mass of volatiles set in MainChamber in real model
    m_h2o_melt = M_h2o / (V * rho_m * (1 - eps_x0))
    m_co2_melt = M_co2 / (V * rho_m * (1 - eps_x0))
    m_eq_max = exsolve_meq(composition, P, T, 0.0)
    if m_h2o_melt > m_eq_max
        phase = 3
    else
        C_co2_sat = exsolve3(composition, P, T, m_h2o_melt)[1]
        if m_co2_melt > C_co2_sat
            phase = 3
        else
            phase = 2
        end
    end
    return (; phase, m_co2_melt)
end


function fun(x::Float64, P::Float64, T::Float64, eps_g0::Float64, eps_x0::Float64, V::Float64, rho_m::Float64, rho_g::Float64, mm_co2::Float64, mm_h2o::Float64, M_co2::Float64)::Float64
    sol = real(mco2_dissolved_sat(x, P, T) * (1 - eps_g0 - eps_x0) * V * rho_m + eps_g0 * rho_g * V * x * mm_co2 / (x * mm_co2 + (1 - x) * mm_h2o) - M_co2)
    return sol
end

function solve_X_co2(eps_g0::Float64, X_co2_prev::Float64, P::Float64, T::Float64, eps_x0::Float64, V::Float64, rho_m::Float64, rho_g::Float64, M_co2::Float64, Tol::Float64)::Float64
    c = ConstantValues()
    mm_co2, mm_h2o = c.mm_co2, c.mm_h2o
    fx = ZeroProblem(x -> fun(x, P, T, eps_g0, eps_x0, V, rho_m, rho_g, mm_co2, mm_h2o, M_co2), X_co2_prev)
    X_co20 = solve(fx; xatol=Tol, atol=Tol, maxiters=100)
    if X_co20 < 0 || X_co20 > 1
        X_co20 = -1.0
        return X_co20
    end
    fraction = ParamICFinder().fraction
    Xmean = (1 - fraction) * X_co2_prev + fraction * X_co20
    X_co20 = Xmean
    return X_co20
end

function get_eps_g(composition::Union{Silicic,Mafic}, eps_g_prev::Float64, X_co20::Float64, P::Float64, T::Float64, eps_x0::Float64, V::Float64, rho_m::Float64, rho_g::Float64, M_h2o::Float64, M_co2::Float64)::NamedTuple{(:eps_g0, :mco2_diss),NTuple{2,Float64}}
    mwater_dissolved = meq_water(composition, X_co20, P, T)
    mco2_diss = mco2_dissolved_sat(X_co20, P, T)
    eps_m = 1 - eps_g_prev - eps_x0
    Num = M_co2 - mco2_diss * eps_m * V * rho_m + M_h2o - mwater_dissolved * eps_m * V * rho_m
    Den = rho_g * V
    fraction = ParamICFinder().fraction
    eps_g0 = (1 - fraction) * eps_g_prev + fraction * Num / Den
    return (; eps_g0, mco2_diss)
end
