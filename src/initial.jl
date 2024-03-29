"""
    eos_g(P::Float64, T::Float64)::EosG{Float64}

parametrization of redlich kwong taken from Huber et al. 2010

## Arguments
- `P`: Pressure (Pa)
- `T`: Temperature (K)

## Returns
An instance of the `EosG` struct, containing the following fields:
- `rho_g`: density of the gas (kg/m³).
- `drho_g_dP`: partial derivative of the density with respect to pressure (kg/m³/Pa).
- `drho_g_dT`: partial derivative of the density with respect to temperature (kg/m³/K).
"""
function eos_g(P::Float64, T::Float64)::EosG{Float64}
    eos_g = EosG(P, T)
    return eos_g
end

"""
    eos_g_rho_g(P::Float64, T::Float64)::Float64

Spetialized version of `eos_g` that computes `rho_g` only.

## Arguments
- `P`: Pressure (Pa)
- `T`: Temperature (K)

## Returns
- `rho_g`: density of the gas in kg/m³.
"""
function eos_g_rho_g(P::Float64, T::Float64)::Float64
    ρ = EosG_RhoG(P, T).rho_g
    return ρ
end

"""
    exsolve(composition::Silicic, P::Float64, T::Float64, X_co2::Float64)::NamedTuple{(:meq, :dmeqdP, :dmeqdT, :dmeqdXco2, :C_co2, :dC_co2dP, :dC_co2dT, :dC_co2dXco2), NTuple{8, Float64}}

This script uses Liu et al. (2006) to calculate the solubility of water

## Arguments
- `composition`: `Silicic()`
- `P`: Pressure (Pa)
- `T`: Temperature (K)
- `X_co2`: mole fraction of CO2 in gas.
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

## Arguments
- `composition`: `Mafic()`
- `P`: Pressure (Pa)
- `T`: Temperature (K)
- `X_co2`: mole fraction of CO2 in gas.
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

Spetialized version of `exsolve` that computes `meq` only.

## Arguments
- `composition`: `Silicic()`
- `P`: Pressure (Pa)
- `T`: Temperature (K)
- `X_co2`: mole fraction of CO2 in gas.

## Returns
- `meq`: amount of water
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

Spetialized version of `exsolve` that computes `meq` only.
## Arguments
- `composition`: `Mafic()`
- `P`: Pressure (Pa)
- `T`: Temperature (K)
- `X_co2`: mole fraction of CO2 in gas.

## Returns
- `meq`: amount of water
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

Takes pressure, temperature, and amount of water to solve for the concentration of CO2 and X_CO2 using a Newton Raphson scheme

## Arguments
- `composition`: `Silicic()`
- `P`: pressure (Pa)
- `T`: temperature (K)
- `m_eq`: amount of water

## Returns
A NamedTuple with the following fields:
- `C_co2`: The concentration of CO2 in the melt.
- `X_co2`: The mole fraction of CO2 in the gas.
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

    # Calculate partial pressures of CO2 and water
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

Takes pressure, temperature, and amount of water to solve for the concentration of CO2 and X_CO2 using a Newton Raphson scheme

## Arguments
- `composition`: `Mafic()`
- `P`: pressure (Pa)
- `T`: temperature (K)
- `m_eq`: amount of water

## Returns
A NamedTuple with the following fields:
- `C_co2`: The concentration of CO2 in the melt.
- `X_co2`: The mole fraction of CO2 in the gas.
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

- This function is used within the `find_liq`, `crystal_fraction`, `crystal_fraction_eps_x` function.

## Arguments
- `composition`: `Silicic()`
- `mH2O`: Weight fration of the H2O in magma.
- `mCO2`: Weight fration of the CO2 in magma.
- `P`: Pressure (Pa)
"""
function parameters_melting_curve(
    composition::Silicic, mH2O::Float64, mCO2::Float64, P::Float64
)::NamedTuple{
    (:a, :dadx, :dady, :dadz, :b, :dbdx, :dbdy, :dbdz, :c, :dcdx, :dcdy, :dcdz),
    NTuple{12,Float64},
}
    x, y, z = mH2O, mCO2, P / 1e6
    a, dadx, dady, dadz = dX_dxdydz(composition, "a", x, y, z)
    b, dbdx, dbdy, dbdz = dX_dxdydz(composition, "b", x, y, z)
    c, dcdx, dcdy, dcdz = dX_dxdydz(composition, "c", x, y, z)
    return (; a, dadx, dady, dadz, b, dbdx, dbdy, dbdz, c, dcdx, dcdy, dcdz)
end

"""
    parameters_melting_curve(composition::Mafic, mH2O::Float64, mCO2::Float64, P::Float64)::NamedTuple{(:a, :dadx, :dady, :dadz, :b, :dbdx, :dbdy, :dbdz), NTuple{8, Float64}}

- This function is used within the `find_liq`, `crystal_fraction`, `crystal_fraction_eps_x` function.

## Arguments
- `composition`: `Mafic()`
- `mH2O`: Weight fration of the H2O in magma.
- `mCO2`: Weight fration of the CO2 in magma.
- `P`: Pressure (Pa)
"""
function parameters_melting_curve(
    composition::Mafic, mH2O::Float64, mCO2::Float64, P::Float64
)::NamedTuple{(:a, :dadx, :dady, :dadz, :b, :dbdx, :dbdy, :dbdz),NTuple{8,Float64}}
    x, y, z = mH2O, mCO2, P / 1e6
    a, dadx, dady, dadz = dX_dxdydz(composition, "a", x, y, z)
    b, dbdx, dbdy, dbdz = dX_dxdydz(composition, "b", x, y, z)
    return (; a, dadx, dady, dadz, b, dbdx, dbdy, dbdz)
end

"""
    find_liq(composition::Silicic, water::Float64, co2::Float64, P::Float64, ini_eps_x::Float64)::Float64

Given the weight fractions of H2O and CO2 in magma, and the pressure, finds the temperature at which the magma will crystallize and reach the desired volume fraction of crystals. The function uses the parameters of the melting curve and an error function to estimate the liquidus temperature.

## Arguments
- `composition`: `Silicic()`
- `water`: Weight fration of the H2O in magma.
- `co2`: Weight fration of the CO2 in magma.
- `P`: Pressure (Pa)
- `ini_eps_x`: The starting volumn fraction of crystal.

## Returns
- `Tl`: the estimated liquidus temperature (K).
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

Given the weight fractions of H2O and CO2 in magma, and the pressure, finds the temperature at which the magma will crystallize and reach the desired volume fraction of crystals. The function uses the parameters of the melting curve and an error function to estimate the liquidus temperature.

## Arguments
- `composition`: `Mafic()`
- `water`: Weight fration of the H2O in magma.
- `co2`: Weight fration of the CO2 in magma.
- `P`: Pressure (Pa)
- `ini_eps_x`: The starting volumn fraction of crystal.

## Returns
- `Tl`: the estimated liquidus temperature (K).
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

Calculate crystal volume fraction and its derivatives in magma at given temperature, pressure, and water and CO2 contents.

## Arguments
- `composition`: `Silicic()`
- `T`: Temperature (K)
- `P`: Pressure (Pa)
- `mH2O`: Weight fration of the H2O in magma.
- `mCO2`: Weight fration of the CO2 in magma.

## Returns
A NamedTuple with the following fields:
- `eps_x`: Crystal volume fraction.
- `deps_x_dP`: Derivative of eps_x with respect to pressure.
- `deps_x_dT`: Derivative of eps_x with respect to temperature.
- `deps_x_deps_g`: Derivative of eps_x with respect to gas volume fraction.
- `deps_x_dmco2_t`: Derivative of eps_x with respect to CO2 content.
- `deps_x_dmh2o_t`: Derivative of eps_x with respect to H2O content.
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

Calculate crystal volume fraction and its derivatives in magma at given temperature, pressure, and water and CO2 contents.

## Arguments
- `composition`: `Mafic()`
- `T`: Temperature (K)
- `P`: Pressure (Pa)
- `mH2O`: Weight fration of the H2O in magma.
- `mCO2`: Weight fration of the CO2 in magma.

## Returns
A NamedTuple with the following fields:
- `eps_x`: Crystal volume fraction.
- `deps_x_dP`: Derivative of eps_x with respect to pressure.
- `deps_x_dT`: Derivative of eps_x with respect to temperature.
- `deps_x_deps_g`: Derivative of eps_x with respect to gas volume fraction.
- `deps_x_dmco2_t`: Derivative of eps_x with respect to CO2 content.
- `deps_x_dmh2o_t`: Derivative of eps_x with respect to H2O content.
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

"""
    crystal_fraction_eps_x(composition::Silicic, T::Float64, P::Float64, mH2O::Float64, mCO2::Float64)::Float64

- Calculates crystal volume fraction at given temperature, pressure, and water and CO2 contents.
- Spetialized version of `crystal_fraction` that computes `eps_x` only.

## Returns
- `eps_x`: Crystal volume fraction.
"""
function crystal_fraction_eps_x(
    composition::Silicic, T::Float64, P::Float64, mH2O::Float64, mCO2::Float64
)::Float64
    T = T - 273.0
    nt = parameters_melting_curve(composition, 100 * mH2O, 100 * mCO2, P)
    eps_x = nt.a * erfc(nt.b * (T - nt.c))
    return eps_x
end

"""
    crystal_fraction_eps_x(composition::Mafic, T::Float64, P::Float64, mH2O::Float64, mCO2::Float64)::Float64

- Calculates crystal volume fraction at given temperature, pressure, and water and CO2 contents.
- Spetialized version of `crystal_fraction`` that computes `eps_x` only.

## Returns
- `eps_x`: Crystal volume fraction.
"""
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
