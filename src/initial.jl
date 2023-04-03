"""
    eos_g(P::Float64, T::Float64)::EosG{Float64}

parametrization of redlich kwong taken from Huber et al. 2010

# Arguments
-`P`: Pressure (Pa)
-`T`: Temperature (K)
"""
function eos_g(P::Float64, T::Float64)::EosG{Float64}
    eos_g = EosG(P, T)
    return eos_g
end

"""
    eos_g_rho_g(P::Float64, T::Float64)::Float64

Spetialized version of eos_g that computes `rho_g` only.
"""
function eos_g_rho_g(P::Float64, T::Float64)::Float64
    ρ = EosG_RhoG(P, T).rho_g
    return ρ
end

"""
    exsolve(composition::Silicic, P::Float64, T::Float64, X_co2::Float64)::NamedTuple{(:meq, :dmeqdP, :dmeqdT, :dmeqdXco2, :C_co2, :dC_co2dP, :dC_co2dT, :dC_co2dXco2), NTuple{8, Float64}}

This script uses Liu et al. (2006) to calculate the solubility of water

# Arguments
-`P`: represents pressure
-`T`: represents the temperature in some units
-`X_co2`: mole fraction of CO2 in gas.
"""
function exsolve(
    composition::Silicic, P::Float64, T::Float64, X_co2::Float64
)::NamedTuple{
    (:meq, :dmeqdP, :dmeqdT, :dmeqdXco2, :C_co2, :dC_co2dP, :dC_co2dT, :dC_co2dXco2),
    NTuple{8,Float64},
}
    # Henry's law
    # partial pressures of CO2 and Water
    P = P / 1e6
    Pc, dPcdP, dPcdXco2 = P * X_co2, X_co2, P
    Pw, dPwdP, dPwdXco2 = P * (1 - X_co2), 1 - X_co2, -P
    meq, dmeqdT, dmeqdP, dmeqdXco2 = build_meq(
        composition, Pw, Pc, T, dPwdP, dPcdP, dPwdXco2, dPcdXco2
    )
    # coefficients for CO2 partitioning
    C_co2, dC_co2dT, dC_co2dP, dC_co2dXco2 = build_co2(
        Pw, Pc, T, dPwdP, dPcdP, dPwdXco2, dPcdXco2
    )
    return (; meq, dmeqdP, dmeqdT, dmeqdXco2, C_co2, dC_co2dP, dC_co2dT, dC_co2dXco2)
end

"""
    exsolve(composition::Mafic, P::Float64, T::Float64, X_co2::Float64)::NamedTuple{(:meq, :dmeqdP, :dmeqdT, :dmeqdXco2, :C_co2, :dC_co2dP, :dC_co2dT, :dC_co2dXco2), NTuple{8, Float64}}

This script uses Liu et al. (2006) to calculate the solubility of water

# Arguments
-`P`: represents pressure
-`T`: represents the temperature in some units
-`X_co2`: mole fraction of CO2 in gas.
"""
function exsolve(
    composition::Mafic, P::Float64, T::Float64, X_co2::Float64
)::NamedTuple{
    (:meq, :dmeqdP, :dmeqdT, :dmeqdXco2, :C_co2, :dC_co2dP, :dC_co2dT, :dC_co2dXco2),
    NTuple{8,Float64},
}
    # Henry's law
    # partial pressures of CO2 and Water
    P = P / 1e6
    Pc, dPcdP, dPcdXco2 = P * X_co2, X_co2, P
    Pw, dPwdP, dPwdXco2 = P * (1 - X_co2), 1 - X_co2, -P
    meq, dmeqdT, dmeqdP, dmeqdXco2 = build_meq(composition, P, T, X_co2)
    # coefficients for CO2 partitioning
    C_co2, dC_co2dT, dC_co2dP, dC_co2dXco2 = build_co2(
        Pw, Pc, T, dPwdP, dPcdP, dPwdXco2, dPcdXco2
    )
    return (; meq, dmeqdP, dmeqdT, dmeqdXco2, C_co2, dC_co2dP, dC_co2dT, dC_co2dXco2)
end

"""
    exsolve_meq(composition::Silicic, P::Float64, T::Float64, X_co2::Float64)::Float64

This script uses Liu et al. (2006) to calculate the solubility of water

# Arguments
-`P`: represents pressure
-`T`: represents the temperature in some units
-`X_co2`: mole fraction of CO2 in gas.
"""
function exsolve_meq(composition::Silicic, P::Float64, T::Float64, X_co2::Float64)::Float64
    P = P / 1e6
    Pc, Pw = P * X_co2, P * (1 - X_co2)
    @unpack b1, b2, b3, b4, b5, b6 = ExsolveSilicic()
    meq = meq_silicic(Pw, Pc, T, b1, b2, b3, b4, b5, b6)
    return meq
end

"""
    exsolve_meq(composition::Mafic, P::Float64, T::Float64, X_co2::Float64)::Float64

This script uses Liu et al. (2006) to calculate the solubility of water

# Arguments
-`P`: represents pressure
-`T`: represents the temperature in some units
-`X_co2`: mole fraction of CO2 in gas.
"""
function exsolve_meq(composition::Mafic, P::Float64, T::Float64, X_co2::Float64)::Float64
    P = P / 1e6
    T_C = T - 273.15
    @unpack b1, b2, b3, b4, b5, b6, b7, b8, b9, b10 = ExsolveMafic()
    meq = meq_mafic(P, T_C, X_co2, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
    return meq
end

"""
    exsolve3(composition::Silicic, P::Float64, T::Float64, m_eq::Float64)::NamedTuple{(:C_co2, :X_co2), NTuple{2, Float64}}

Takes pressure, temperature, and amount of water to solve for the concentration of CO2 and X_CO2 (basically, goes the other direction compared to exsolve) using a Newton Raphson scheme

# Arguments
-`P`: pressure (Pa)
-`T`: temperature (K)
-`m_eq`: amount of water
"""
function exsolve3(
    composition::Silicic, P::Float64, T::Float64, m_eq::Float64
)::NamedTuple{(:C_co2, :X_co2),NTuple{2,Float64}}
    # convert to MPa and Celsius
    P = P / 1e6
    m_eq = m_eq * 1e2

    f(Xc_prev::Float64)::Float64 = water(composition, P, T, Xc_prev, m_eq)
    f_prime(Xc_prev::Float64)::Float64 = dwater_dx(composition, P, T, Xc_prev)

    X_co2 = solve_NR(f, f_prime, 1e-10, 1e2, 1e-2)

    # partial pressures of CO2 and Water
    Pc = P * X_co2
    Pw = P * (1 - X_co2)

    # function & coefficients from Liu et al 2005
    @unpack c1, c2, c3, c4 = Co2PartitionCoeff()
    C_co2 = Pc * (c1 + c2 * Pw) / T + Pc * (c3 * Complex(Pw)^0.5 + c4 * Complex(Pw)^1.5)
    C_co2 = real(C_co2 * 1e-6)
    return (; C_co2, X_co2)
end

"""
    exsolve3(composition::Mafic, P::Float64, T::Float64, m_eq::Float64)::NamedTuple{(:C_co2, :X_co2), NTuple{2, Float64}}

Takes pressure, temperature, and amount of water to solve for the concentration of CO2 and X_CO2 (basically, goes the other direction compared to exsolve.m) using a Newton Raphson scheme

# Arguments
-`P`: pressure (Pa)
-`T`: temperature (K)
-`m_eq`: amount of water
"""
function exsolve3(
    composition::Mafic, P::Float64, T::Float64, m_eq::Float64
)::NamedTuple{(:C_co2, :X_co2),NTuple{2,Float64}}
    # convert to MPa and Celsius
    P = P / 1e6
    T = T - 273.15
    m_eq = m_eq * 1e2

    f(Xc_prev::Float64)::Float64 = water(composition, P, T, Xc_prev, m_eq)
    f_prime(Xc_prev::Float64)::Float64 = dwater_dx(composition, P, T, Xc_prev)

    X_co2 = solve_NR(f, f_prime, 1e-10, 1e2, 1e-2)

    # partial pressures of CO2 and Water
    Pc = P * X_co2
    Pw = P * (1 - X_co2)
    T = T + 273.15  # convert back to Kelvin because that's what the Liu needs

    # function & coefficients from Liu et al 2005
    @unpack c1, c2, c3, c4 = Co2PartitionCoeff()
    C_co2 = Pc * (c1 + c2 * Pw) / T + Pc * (c3 * Complex(Pw)^0.5 + c4 * Complex(Pw)^1.5)
    C_co2 = real(C_co2 * 1e-6)
    return (; C_co2, X_co2)
end

"""
    parameters_melting_curve(composition::Silicic, mH2O::Float64, mCO2::Float64, P::Float64)::NamedTuple{(:a, :dadx, :dady, :dadz, :b, :dbdx, :dbdy, :dbdz, :c, :dcdx, :dcdy, :dcdz), NTuple{12, Float64}}

# Arguments
-`mH2O`: Weight fration of the H2O in magma.
-`mCO2`: Weight fration of the CO2 in magma.
-`P`: Pressure (Pa)
"""
function parameters_melting_curve(
    composition::Silicic, mH2O::Float64, mCO2::Float64, P::Float64
)::NamedTuple{
    (:a, :dadx, :dady, :dadz, :b, :dbdx, :dbdy, :dbdz, :c, :dcdx, :dcdy, :dcdz),
    NTuple{12,Float64},
}
    x = mH2O
    y = mCO2
    z = P / 1e6

    a =
        0.36 - 0.02 * x - 0.06 * y +
        8.6E-04 * z +
        0.0024 * x * y +
        6.27E-05 * x * z +
        3.57E-05 * y * z - 0.0026 * x^2 + 0.003 * y^2 - 1.16E-06 * z^2
    dadx = 100 * (-0.02 + 0.0024 * y + 6.27E-05 * z - 2 * 0.0026 * x)
    dady = 100 * (-0.06 + 0.0024 * x + 3.57E-05 * z + 2 * 0.003 * y)
    dadz = 1e-6 * (8.6E-04 + 6.27E-05 * x + 3.57E-05 * y - 2 * 1.16E-06 * z)

    b =
        0.0071 + 0.0049 * x + 0.0043 * y - 4.08E-05 * z - 7.85E-04 * x * y -
        1.3E-05 * x * z +
        3.97E-06 * y * z +
        6.29E-04 * x^2 - 0.0025 * y^2 + 8.51E-08 * z^2
    dbdx = 100 * (0.0049 - 7.85E-04 * y - 1.3E-05 * z + 2 * 6.29E-04 * x)
    dbdy = 100 * (0.0043 - 7.85E-04 * x + 3.97E-06 * z - 2 * 0.0025 * y)
    dbdz = 1e-6 * (-4.08E-05 - 1.3E-05 * x + 3.97E-06 * y + 2 * 8.51E-08 * z)

    c =
        863.09 - 36.9 * x + 48.81 * y - 0.17 * z - 1.52 * x * y - 0.04 * x * z -
        0.04 * y * z + 4.57 * x^2 - 7.79 * y^2 + 4.65E-04 * z^2
    dcdx = 100 * (-36.9 - 1.52 * y - 0.04 * z + 2 * 4.57 * x)
    dcdy = 100 * (48.81 - 1.52 * x - 0.04 * z - 2 * 7.79 * y)
    dcdz = 1e-6 * (-0.17 - 0.04 * x - 0.04 * y + 2 * 4.65E-04 * z)
    return (; a, dadx, dady, dadz, b, dbdx, dbdy, dbdz, c, dcdx, dcdy, dcdz)
end

"""
    parameters_melting_curve(composition::Mafic, mH2O::Float64, mCO2::Float64, P::Float64)::NamedTuple{(:a, :dadx, :dady, :dadz, :b, :dbdx, :dbdy, :dbdz), NTuple{8, Float64}}

# Arguments
-`mH2O`: Weight fration of the H2O in magma.
-`mCO2`: Weight fration of the CO2 in magma.
-`P`: Pressure (Pa)
"""
function parameters_melting_curve(
    composition::Mafic, mH2O::Float64, mCO2::Float64, P::Float64
)::NamedTuple{(:a, :dadx, :dady, :dadz, :b, :dbdx, :dbdy, :dbdz),NTuple{8,Float64}}
    x = mH2O
    y = mCO2
    z = P / 1e6

    intercept = -0.0106007180771044
    H2Ocoeff = 0.00642879427079997
    CO2coeff = -0.000362698886591479
    Pcoeff = 6.33762356329763e-06
    H2OxCO2coeff = -0.0000409319695736593
    H2OxPcoeff = -0.0000020971242285322
    CO2xPcoeff = 3.66354014084072e-07
    H2Osquarecoeff = -0.00127225661031757
    CO2squarecoeff = 0.000219802992993448
    Psquarecoeff = -1.4241625041626E-09

    a =
        intercept +
        H2Ocoeff * x +
        CO2coeff * y +
        Pcoeff * z +
        H2OxCO2coeff * x * y +
        H2OxPcoeff * x * z +
        CO2xPcoeff * y * z +
        H2Osquarecoeff * x^2 +
        CO2squarecoeff * y^2 +
        Psquarecoeff * z^2
    dadx = 100 * (H2Ocoeff + H2OxCO2coeff * y + H2OxPcoeff * z + 2 * H2Osquarecoeff * x)
    dady = 100 * (CO2coeff + H2OxCO2coeff * x + CO2xPcoeff * z + 2 * CO2squarecoeff * y)
    dadz = 1e-6 * (Pcoeff + H2OxPcoeff * x + CO2xPcoeff * y + 2 * Psquarecoeff * z)

    # b value
    intercept = 12.1982401917454
    H2Ocoeff = -7.49690626527448
    CO2coeff = 0.398381500262876
    Pcoeff = -0.00632911929609247
    H2OxCO2coeff = 0.0571369994114008
    H2OxPcoeff = 0.00216190962922558
    CO2xPcoeff = -0.000409810092770206
    H2Osquarecoeff = 1.48907741502382
    CO2squarecoeff = -0.251451720536687
    Psquarecoeff = 1.36630369630388e-06

    b =
        intercept +
        H2Ocoeff * x +
        CO2coeff * y +
        Pcoeff * z +
        H2OxCO2coeff * x * y +
        H2OxPcoeff * x * z +
        CO2xPcoeff * y * z +
        H2Osquarecoeff * x^2 +
        CO2squarecoeff * y^2 +
        Psquarecoeff * z^2
    dbdx = 100 * (H2Ocoeff + H2OxCO2coeff * y + H2OxPcoeff * z + 2 * H2Osquarecoeff * x)
    dbdy = 100 * (CO2coeff + H2OxCO2coeff * x + CO2xPcoeff * z + 2 * CO2squarecoeff * y)
    dbdz = 1e-6 * (Pcoeff + H2OxPcoeff * x + CO2xPcoeff * y + 2 * Psquarecoeff * z)
    return (; a, dadx, dady, dadz, b, dbdx, dbdy, dbdz)
end

"""
    find_liq(composition::Silicic, water::Float64, co2::Float64, P::Float64, ini_eps_x::Float64)::Float64

# Arguments
-`water`: Weight fration of the H2O in magma.
-`co2`: Weight fration of the CO2 in magma.
-`P`: Pressure (Pa)
-`ini_eps_x`: The starting volumn fraction of crystal.
"""
function find_liq(
    composition::Silicic, water::Float64, co2::Float64, P::Float64, ini_eps_x::Float64
)::Float64
    nt = parameters_melting_curve(composition, 100 * water, 100 * co2, P)
    f(x) = nt.a * erfc(nt.b * (x - nt.c)) - ini_eps_x
    x0 = 1000
    Tl = fzero(f, (0, x0); maxevals=100)
    return Tl + 273.15
end

"""
    find_liq(composition::Mafic, water::Float64, co2::Float64, P::Float64, ini_eps_x::Float64)::Float64

# Arguments
-`water`: Weight fration of the H2O in magma.
-`co2`: Weight fration of the CO2 in magma.
-`P`: Pressure (Pa)
-`ini_eps_x`: The starting volumn fraction of crystal.
"""
function find_liq(
    composition::Mafic, water::Float64, co2::Float64, P::Float64, ini_eps_x::Float64
)::Float64
    nt = parameters_melting_curve(composition, 100 * water, 100 * co2, P)
    Tl = (ini_eps_x - nt.b) / nt.a
    return Tl + 273.15
end

"""
    crystal_fraction(composition::Silicic, T::Float64, P::Float64, mH2O::Float64, mCO2::Float64)::NamedTuple{(:eps_x, :deps_x_dP, :deps_x_dT, :deps_x_deps_g, :deps_x_dmco2_t, :deps_x_dmh2o_t), NTuple{6, Float64}}

# Arguments
-`T`: Temperature (K)
-`P`: Pressure (Pa)
-`mH2O`: Weight fration of the H2O in magma.
-`mCO2`: Weight fration of the CO2 in magma.
"""
function crystal_fraction(
    composition::Silicic, T::Float64, P::Float64, mH2O::Float64, mCO2::Float64
)::NamedTuple{
    (:eps_x, :deps_x_dP, :deps_x_dT, :deps_x_deps_g, :deps_x_dmco2_t, :deps_x_dmh2o_t),
    NTuple{6,Float64},
}
    # NEW VERSION WITH SAGE's PARAMETERIZATION
    T = T - 273
    a, dadx, dady, dadz, b, dbdx, dbdy, dbdz, c, dcdx, dcdy, dcdz = parameters_melting_curve(
        composition, 100 * mH2O, 100 * mCO2, P
    )
    eps_x = a * erfc(b * (T - c))
    # derivatives
    deps_x_deps_g = -1.0
    deps_x_dT = -2 * a * b * exp(-b^2 * (T - c)^2) / sqrt(pi)
    deps_x_dP =
        dadz * erfc(b * (T - c)) -
        2 * a * (T - c) / sqrt(pi) * exp(-b^2 * (T - c)^2) * dbdz +
        2 * a * b / sqrt(pi) * exp(-b^2 * (T - c)^2) * dcdz
    deps_x_dmco2_t =
        dady * erfc(b * (T - c)) -
        2 * a * (T - c) / sqrt(pi) * exp(-b^2 * (T - c)^2) * dbdy +
        2 * a * b / sqrt(pi) * exp(-b^2 * (T - c)^2) * dcdy
    deps_x_dmh2o_t =
        dadx * erfc(b * (T - c)) -
        2 * a * (T - c) / sqrt(pi) * exp(-b^2 * (T - c)^2) * dbdx +
        2 * a * b / sqrt(pi) * exp(-b^2 * (T - c)^2) * dcdx
    return (; eps_x, deps_x_dP, deps_x_dT, deps_x_deps_g, deps_x_dmco2_t, deps_x_dmh2o_t)
end

"""
    crystal_fraction(composition::Mafic, T::Float64, P::Float64, mH2O::Float64, mCO2::Float64)::NamedTuple{(:eps_x, :deps_x_dP, :deps_x_dT, :deps_x_deps_g, :deps_x_dmco2_t, :deps_x_dmh2o_t), NTuple{6, Float64}}

# Arguments
-`T`: Temperature (K)
-`P`: Pressure (Pa)
-`mH2O`: Weight fration of the H2O in magma.
-`mCO2`: Weight fration of the CO2 in magma.
"""
function crystal_fraction(
    composition::Mafic, T::Float64, P::Float64, mH2O::Float64, mCO2::Float64
)::NamedTuple{
    (:eps_x, :deps_x_dP, :deps_x_dT, :deps_x_deps_g, :deps_x_dmco2_t, :deps_x_dmh2o_t),
    NTuple{6,Float64},
}
    # NEW VERSION WITH SAGE's PARAMETERIZATION
    T = T - 273
    a, dadx, dady, dadz, b, dbdx, dbdy, dbdz = parameters_melting_curve(
        composition, 100 * mH2O, 100 * mCO2, P
    )
    eps_x = a * T + b
    # derivatives
    deps_x_deps_g = -1.0
    deps_x_dT = a

    deps_x_dmco2_t = dady * T + dbdy
    deps_x_dmh2o_t = dadx * T + dbdx
    deps_x_dP = dadz * T + dbdz

    if eps_x < 0 || eps_x > 1
        eps_x = 0.0
        deps_x_deps_g = 0.0
        deps_x_dT = 0.0
    end
    return (; eps_x, deps_x_dP, deps_x_dT, deps_x_deps_g, deps_x_dmco2_t, deps_x_dmh2o_t)
end

function crystal_fraction_eps_x(
    composition::Silicic, T::Float64, P::Float64, mH2O::Float64, mCO2::Float64
)::Float64
    T = T - 273.0
    nt = parameters_melting_curve(composition, 100 * mH2O, 100 * mCO2, P)
    eps_x = nt.a * erfc(nt.b * (T - nt.c))
    return eps_x
end

function crystal_fraction_eps_x(
    composition::Mafic, T::Float64, P::Float64, mH2O::Float64, mCO2::Float64
)::Float64
    T = T - 273.0
    nt = parameters_melting_curve(composition, 100 * mH2O, 100 * mCO2, P)
    eps_x = nt.a * T + nt.b
    if eps_x < 0 || eps_x > 1
        eps_x = 0.0
    end
    return eps_x
end
