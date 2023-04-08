# For function exsolve
struct ExsolveSilicic{T}
    b1::T
    b2::T
    b3::T
    b4::T
    b5::T
    b6::T
    ExsolveSilicic() = new{Float64}(354.94, 9.623, -1.5223, 1.2439e-3, -1.084e-4, -1.362e-5)
end

struct ExsolveMafic{T}
    b1::T
    b2::T
    b3::T
    b4::T
    b5::T
    b6::T
    b7::T
    b8::T
    b9::T
    b10::T
    function ExsolveMafic()
        return new{Float64}(
            2.99622526644026,
            0.00322422830627781,
            -9.1389095360385,
            0.0336065247530767,
            0.00747236662935722,
            -0.0000150329805347769,
            -0.01233608521548,
            -4.14842647942619e-6,
            -0.655454303068124,
            -7.35270395041104e-6,
        )
    end
end

struct Co2PartitionCoeff{T}
    c1::T
    c2::T
    c3::T
    c4::T
    Co2PartitionCoeff() = new{Float64}(5668, -55.99, 0.4133, 2.041e-3)
end

function meq_silicic(
    Pw::T, Pc::T, Temp::T, b1::T, b2::T, b3::T, b4::T, b5::T, b6::T
)::T where {T<:Float64}
    return 1e-2 * real(
        (b1 * Complex(Pw)^0.5 + b2 * Pw + b3 * Complex(Pw)^1.5) / Temp +
        b4 * Complex(Pw)^1.5 +
        Pc * (b5 * Complex(Pw)^0.5 + b6 * Pw),
    )
end

function dmeqdT_silicic(Pw::T, Temp::T, b1::T, b2::T, b3::T)::T where {T<:Float64}
    return 1e-2 *
           real(-1 * (b1 * Complex(Pw)^0.5 + b2 * Pw + b3 * Complex(Pw)^1.5) * Temp^(-2))
end

function dmeqdP_silicic(
    Pw::T, dPwdP::T, Pc::T, dPcdP::T, Temp::T, b1::T, b2::T, b3::T, b4::T, b5::T, b6::T
)::T where {T<:Float64}
    return 1e-8 * real(
        dPwdP * (
            (1 / Temp) * (0.5 * b1 * Complex(Pw)^(-0.5) + b2 + 1.5 * b3 * Complex(Pw)^0.5) +
            1.5 * b4 * Complex(Pw)^0.5 +
            Pc * (0.5 * b5 * Complex(Pw)^(-0.5) + b6)
        ) + dPcdP * (b5 * Complex(Pw)^0.5 + b6 * Pw),
    )
end

function dmeqdXco2_silicic(
    Pw::T,
    dPwdXco2::T,
    Pc::T,
    dPcdXco2::T,
    Temp::T,
    b1::T,
    b2::T,
    b3::T,
    b4::T,
    b5::T,
    b6::T,
)::T where {T<:Float64}
    return 1e-2 * real(
        dPwdXco2 * (
            (1 / Temp) * (0.5 * b1 * Complex(Pw)^(-0.5) + b2 + 1.5 * b3 * Complex(Pw)^0.5) +
            1.5 * b4 * Complex(Pw)^0.5 +
            Pc * (0.5 * b5 * Complex(Pw)^(-0.5) + b6)
        ) + dPcdXco2 * (b5 * Complex(Pw)^0.5 + b6 * Pw),
    )
end

function meq_mafic(
    P::T,
    T_C::T,
    X_co2::T,
    b1::T,
    b2::T,
    b3::T,
    b4::T,
    b5::T,
    b6::T,
    b7::T,
    b8::T,
    b9::T,
    b10::T,
)::T where {T<:Float64}
    return 1e-2 * real(
        b1 +
        b2 * T_C +
        b3 * X_co2 +
        b4 * P +
        b5 * T_C * X_co2 +
        b6 * T_C * P +
        b7 * X_co2 * P +
        b8 * T_C^2 +
        b9 * Complex(X_co2)^2 +
        b10 * Complex(P)^2,
    )
end

function dmeqdT_mafic(
    P::T, T_C::T, X_co2::T, b2::T, b5::T, b6::T, b8::T
)::T where {T<:Float64}
    return 1e-2 * real(b2 + b5 * X_co2 + b6 * P + 2 * b8 * T_C)
end

function dmeqdP_mafic(
    P::T, T_C::T, X_co2::T, b4::T, b6::T, b7::T, b10::T
)::T where {T<:Float64}
    return 1e-8 * real(b4 + b6 * T_C + b7 * X_co2 + 2 * b10 * P)
end

function dmeqdXco2_mafic(
    P::T, T_C::T, X_co2::T, b3::T, b5::T, b7::T, b9::T
)::T where {T<:Float64}
    return 1e-2 * real(b3 + b5 * T_C + b7 * P + 2 * b9 * X_co2)
end

function C_co2_f(Pw::T, Pc::T, Temp::T, c1::T, c2::T, c3::T, c4::T)::T where {T<:Float64}
    return 1e-6 * real(
        Pc * (c1 + c2 * Pw) / Temp + Pc * (c3 * Complex(Pw)^0.5 + c4 * Complex(Pw)^1.5)
    )
end

function dC_co2dT_f(Pw::T, Pc::T, Temp::T, c1::T, c2::T)::T where {T<:Float64}
    return 1e-6 * real(-1 * Pc * (c1 + c2 * Pw) * Complex(Temp)^(-2))
end

function dC_co2dP_f(
    Pw::T, Pc::T, Temp::T, dPwdP::T, dPcdP::T, c1::T, c2::T, c3::T, c4::T
)::T where {T<:Float64}
    return 1e-12 * real(
        Complex(Temp)^(-1) * (dPcdP * (c1 + c2 * Pw) + Pc * (c2 * dPwdP)) +
        dPcdP * (c3 * Complex(Pw)^0.5 + c4 * Complex(Pw)^1.5) +
        Pc * (0.5 * c3 * Complex(Pw)^(-0.5) * dPwdP + 1.5 * c4 * Complex(Pw)^0.5 * dPwdP),
    )
end

function dC_co2dXco2_f(
    Pw::T, Pc::T, Temp::T, dPwdXco2::T, dPcdXco2::T, c1::T, c2::T, c3::T, c4::T
)::T where {T<:Float64}
    return 1e-6 * real(
        Complex(Temp)^(-1) * (dPcdXco2 * (c1 + c2 * Pw) + Pc * (c2 * dPwdXco2)) +
        dPcdXco2 * (c3 * Complex(Pw)^0.5 + c4 * Complex(Pw)^1.5) +
        Pc *
        (0.5 * c3 * Complex(Pw)^(-0.5) * dPwdXco2 + 1.5 * c4 * Complex(Pw)^0.5 * dPwdXco2),
    )
end

# For function exsolve
"""
    build_meq(composition::Silicic, Pw::T, Pc::T, Temp::T, dPwdP::T, dPcdP::T, dPwdXco2::T, dPcdXco2::T)::NamedTuple{(:meq, :dmeqdT, :dmeqdP, :dmeqdXco2),NTuple{4,T}} where {T<:Float64}

- This function is used within the `exsolve` function.
"""
function build_meq(
    composition::Silicic,
    Pw::T,
    Pc::T,
    Temp::T,
    dPwdP::T,
    dPcdP::T,
    dPwdXco2::T,
    dPcdXco2::T,
)::NamedTuple{(:meq, :dmeqdT, :dmeqdP, :dmeqdXco2),NTuple{4,T}} where {T<:Float64}
    @unpack b1, b2, b3, b4, b5, b6 = ExsolveSilicic()
    meq = meq_silicic(Pw, Pc, Temp, b1, b2, b3, b4, b5, b6)
    dmeqdT = dmeqdT_silicic(Pw, Temp, b1, b2, b3)
    dmeqdP = dmeqdP_silicic(Pw, dPwdP, Pc, dPcdP, Temp, b1, b2, b3, b4, b5, b6)
    dmeqdXco2 = dmeqdXco2_silicic(Pw, dPwdXco2, Pc, dPcdXco2, Temp, b1, b2, b3, b4, b5, b6)
    return (; meq, dmeqdT, dmeqdP, dmeqdXco2)
end

"""
    build_meq(composition::Mafic, P::T, Temp::T, X_co2::T)::NamedTuple{(:meq, :dmeqdT, :dmeqdP, :dmeqdXco2),NTuple{4,T}} where {T<:Float64}

- This function is used within the `exsolve` function.
"""
function build_meq(
    composition::Mafic, P::T, Temp::T, X_co2::T
)::NamedTuple{(:meq, :dmeqdT, :dmeqdP, :dmeqdXco2),NTuple{4,T}} where {T<:Float64}
    T_C = Temp - 273.15
    @unpack b1, b2, b3, b4, b5, b6, b7, b8, b9, b10 = ExsolveMafic()
    meq = meq_mafic(P, T_C, X_co2, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
    dmeqdT = dmeqdT_mafic(P, T_C, X_co2, b2, b5, b6, b8)
    dmeqdP = dmeqdP_mafic(P, T_C, X_co2, b4, b6, b7, b10)
    dmeqdXco2 = dmeqdXco2_mafic(P, T_C, X_co2, b3, b5, b7, b9)
    return (; meq, dmeqdT, dmeqdP, dmeqdXco2)
end

"""
    build_co2(Pw::T, Pc::T, Temp::T, dPwdP::T, dPcdP::T, dPwdXco2::T, dPcdXco2::T)::NamedTuple{(:C_co2, :dC_co2dT, :dC_co2dP, :dC_co2dXco2),NTuple{4,T}} where {T<:Float64}

- This function is used within the `exsolve` function.
"""
function build_co2(
    Pw::T, Pc::T, Temp::T, dPwdP::T, dPcdP::T, dPwdXco2::T, dPcdXco2::T
)::NamedTuple{(:C_co2, :dC_co2dT, :dC_co2dP, :dC_co2dXco2),NTuple{4,T}} where {T<:Float64}
    @unpack c1, c2, c3, c4 = Co2PartitionCoeff()
    C_co2 = C_co2_f(Pw, Pc, Temp, c1, c2, c3, c4)
    dC_co2dT = dC_co2dT_f(Pw, Pc, Temp, c1, c2)
    dC_co2dP = dC_co2dP_f(Pw, Pc, Temp, dPwdP, dPcdP, c1, c2, c3, c4)
    dC_co2dXco2 = dC_co2dXco2_f(Pw, Pc, Temp, dPwdXco2, dPcdXco2, c1, c2, c3, c4)
    return (; C_co2, dC_co2dT, dC_co2dP, dC_co2dXco2)
end

# For function eos_g
rho_g(P::Float64, T::Float64)::Float64 = real(
    1e3 * (
        -112.528 * Complex(T - 273.15)^-0.381 +
        127.811 * Complex(P * 1e-5)^-1.135 +
        112.04 * Complex(T - 273.15)^-0.411 * Complex(P * 1e-5)^0.033
    ),
)
drho_g_dP(P::Float64, T::Float64)::Float64 = real(
    1e-2 * (
        (-1.135) * 127.811 * Complex(P * 1e-5)^-2.135 +
        0.033 * 112.04 * Complex(T - 273.15)^-0.411 * Complex(P * 1e-5)^-0.967
    ),
)
drho_g_dT(P::Float64, T::Float64)::Float64 = real(
    1e3 * (
        (-0.381) * (-112.528) * Complex(T - 273.15)^-1.381 +
        (-0.411) * 112.04 * Complex(T - 273.15)^-1.411 * Complex(P * 1e-5)^0.033
    ),
)

struct EosG{T}
    rho_g::T
    drho_g_dP::T
    drho_g_dT::T
end
EosG(P, T) = EosG{typeof(P)}(rho_g(P, T), drho_g_dP(P, T), drho_g_dT(P, T))

struct EosG_RhoG{T}
    rho_g::T
end
EosG_RhoG(P, T) = EosG_RhoG{typeof(P)}(rho_g(P, T))

# For function exsolve3
@with_kw struct Exsolve3Silicic{T}
    h1::T = 354.94
    h2::T = 9.623
    h3::T = -1.5223
    h4::T = 0.0012439
    h5::T = -1.084e-4
    h6::T = -1.362e-5
end

@with_kw struct Exsolve3Mafic{T}
    h1::T = 2.99622526644026
    h2::T = 0.00322422830627781
    h3::T = -9.1389095360385
    h4::T = 0.0336065247530767
    h5::T = 0.00747236662935722
    h6::T = -0.0000150329805347769
    h7::T = -0.01233608521548
    h8::T = -4.14842647942619e-6
    h9::T = -0.655454303068124
    h10::T = -7.35270395041104e-6
end

# Water Paritioning Function
"""
    water(composition::Silicic, p::T, t::T, x::T, c::T)::T where {T<:Float64}

- Water Paritioning Function
- This function is used within the `exsolve3` function.

# Arguments:
- `p`: pressure (Pa)
- `t`: temperature (K)
- `x`: The previous mole fraction of CO2 (X_CO2)
- `c`: amount of water
"""
function water(composition::Silicic, p::T, t::T, x::T, c::T)::T where {T<:Float64}
    @unpack h1, h2, h3, h4, h5, h6 = Exsolve3Silicic()
    return real(
        (
            h1 * Complex(p * (1 - x))^0.5 +
            h2 * (p * (1 - x)) +
            h3 * Complex(p * (1 - x))^1.5
        ) / t +
        h4 * Complex(p * (1 - x))^1.5 +
        p * x * (h5 * Complex(p * (1 - x))^0.5 + h6 * (p * (1 - x))) - c,
    )
end

"""
    water(composition::Silicic, p::T, t::T, x::T, c::T)::T where {T<:Float64}

- Water Paritioning Function
- This function is used within the `exsolve3` function.

# Arguments:
- `p`: pressure (Pa)
- `t`: temperature (K)
- `x`: The previous mole fraction of CO2 (X_CO2)
- `c`: amount of water
"""
function water(composition::Mafic, p::T, t::T, x::T, c::T)::T where {T<:Float64}
    @unpack h1, h2, h3, h4, h5, h6, h7, h8, h9, h10 = Exsolve3Mafic()
    return real(
        h1 +
        h2 * t +
        h3 * x +
        h4 * p +
        h5 * t * x +
        h6 * t * p +
        h7 * x * p +
        h8 * Complex(t)^2 +
        h9 * x^2 +
        h10 * Complex(p)^2 - c,
    )
end

# Derivative of Water wrt Xco2
"""
    dwater_dx(composition::Silicic, p::T, t::T, x::T)::T where {T<:Float64}

- Function for the derivative of water with reference to X_CO2
- This function is used within the `exsolve3` function.

# Arguments:
- `p`: pressure (Pa)
- `t`: temperature (K)
- `x`: The previous mole fraction of CO2 (X_CO2)
"""
function dwater_dx(composition::Silicic, p::T, t::T, x::T)::T where {T<:Float64}
    @unpack h1, h2, h3, h4, h5, h6 = Exsolve3Silicic()
    return real(
        -p * (
            (1 / t) * (
                0.5 * h1 * Complex(p * (1 - x))^(-0.5) +
                h2 +
                1.5 * h3 * Complex(p * (1 - x))^0.5
            ) +
            1.5 * h4 * Complex(p * (1 - x))^0.5 +
            (p * x) * (0.5 * h5 * Complex(p * (1 - x))^(-0.5) + h6)
        ) + p * (h5 * Complex(p * (1 - x))^0.5 + h6 * (p * (1 - x))),
    )
end

"""
    dwater_dx(composition::Mafic, p::T, t::T, x::T)::T where {T<:Float64}

- Derivative of Water Partitioning Function with respect to X_CO2
- This function is used within the `exsolve3` function.

# Arguments:
- `p`: pressure (Pa)
- `t`: temperature (K)
- `x`: The previous mole fraction of CO2 (X_CO2)
"""
function dwater_dx(composition::Mafic, p::T, t::T, x::T)::T where {T<:Float64}
    @unpack h3, h5, h7, h9 = Exsolve3Mafic()
    return h3 + h5 * t + h7 * p + 2 * h9 * x
end

"""
    solve_NR(f, f_prime, errorTol::T, count_max::T, Xc_initial::T)::T where {T<:Float64}

- Finding mole fraction of CO2 in gas (X_CO2).
- This function is used within the `exsolve3` function.

# Arguments
- `f`: Water Paritioning Function
- `f_prime`: Function for the derivative of water with reference to X_CO2
- `errorTol`: error tolerance, the default value is 1e-10
- `count_max`: Maximum loop count, the default value is 1e2
- `Xc_initial`: initial guess of X_CO2, the dedault value is 1e-2

# Returns
- `X_CO2`: mole fraction of CO2 in gas
"""
function solve_NR(
    f, f_prime, errorTol::T, count_max::T, Xc_initial::T
)::T where {T<:Float64}
    # P,T and inital guesses/values
    Xc_guess = Xc_initial
    Xc_prev = 0.0
    count = 0
    while abs(Xc_prev - Xc_guess) > errorTol && count < count_max
        count = count + 1
        Xc_prev = Xc_guess
        Xc_guess = Xc_prev - (f(Xc_prev) / f_prime(Xc_prev))
    end
    if abs(Xc_prev - Xc_guess) > errorTol && count >= count_max
        Xc_initial = 1e-4
        Xc_guess = Xc_initial
        Xc_prev = 0.0
        count = 0
        while abs(Xc_prev - Xc_guess) > errorTol && count < count_max
            count = count + 1
            Xc_prev = Xc_guess
            Xc_guess = Xc_prev - (f(Xc_prev) / f_prime(Xc_prev))
        end
    end
    while ~isreal(Xc_guess) && Xc_initial <= 1
        Xc_initial = Xc_initial + 0.01
        Xc_prev = 0.0
        while abs(Xc_prev - Xc_guess) > errorTol
            count = count + 1
            Xc_prev = Xc_guess
            Xc_guess = Xc_prev - (f(Xc_prev) / f_prime(Xc_prev))
        end
    end
    X_co2 = Xc_guess
    # MT adding this b/c some super CO2-rich cases make Xc_guess greater than 1
    if X_co2 > 1
        X_co2 = 1.0
    end
    return X_co2
end

# For parameters_melting_curve(Silicic())
@with_kw struct MeltingCurveSilicicA{T}
    c1::T = 0.36
    c2::T = -0.02
    c3::T = -0.06
    c4::T = 8.6e-4
    c5::T = 0.0024
    c6::T = 6.27e-5
    c7::T = 3.57e-5
    c8::T = -0.0026
    c9::T = 0.003
    c10::T = -1.16e-6
end

@with_kw struct MeltingCurveSilicicB{T}
    c1::T = 0.0071
    c2::T = 0.0049
    c3::T = 0.0043
    c4::T = -4.08e-5
    c5::T = -7.85e-4
    c6::T = -1.3e-5
    c7::T = 3.97e-6
    c8::T = 6.29e-4
    c9::T = -0.0025
    c10::T = 8.51e-8
end

@with_kw struct MeltingCurveSilicicC{T}
    c1::T = 863.09
    c2::T = -36.9
    c3::T = 48.81
    c4::T = -0.17
    c5::T = -1.52
    c6::T = -0.04
    c7::T = -0.04
    c8::T = 4.57
    c9::T = -7.79
    c10::T = 4.65e-4
end

# For parameters_melting_curve(Mafic())
@with_kw struct MeltingCurveMaficA{T}
    c1::T = -0.0106007180771044    # intercept
    c2::T = 0.00642879427079997    # H2Ocoeff
    c3::T = -0.000362698886591479  # CO2coeff
    c4::T = 6.33762356329763e-06   # Pcoeff
    c5::T = -0.0000409319695736593 # H2OxCO2coeff
    c6::T = -0.0000020971242285322 # H2OxPcoeff
    c7::T = 3.66354014084072e-07   # CO2xPcoeff
    c8::T = -0.00127225661031757   # H2Osquarecoeff
    c9::T = 0.000219802992993448   # CO2squarecoeff
    c10::T = -1.4241625041626E-09  # Psquarecoeff
end

@with_kw struct MeltingCurveMaficB{T}
    c1::T = 12.1982401917454      # intercept
    c2::T = -7.49690626527448     # H2Ocoeff
    c3::T = 0.398381500262876     # CO2coeff
    c4::T = -0.00632911929609247  # Pcoeff
    c5::T = 0.0571369994114008    # H2OxCO2coeff
    c6::T = 0.00216190962922558   # H2OxPcoeff
    c7::T = -0.000409810092770206 # CO2xPcoeff
    c8::T = 1.48907741502382      # H2Osquarecoeff
    c9::T = -0.251451720536687    # CO2squarecoeff
    c10::T = 1.36630369630388e-06 # Psquarecoeff
end

meltingcurve_dict = Dict(
    Silicic() => Dict(
        "a" => MeltingCurveSilicicA(),
        "b" => MeltingCurveSilicicB(),
        "c" => MeltingCurveSilicicC(),
    ),
    Mafic() => Dict("a" => MeltingCurveMaficA(), "b" => MeltingCurveMaficB()),
)

#=
TODO:
- use more meaningful variable names instead of "a", "b", "c".
=#

"""
    dX_dxdydz(composition::Union{Silicic,Mafic}, s::String, x::T, y::T, z::T)::NamedTuple{(:X, :dXdx, :dXdy, :dXdz),NTuple{4,T}} where {T<:Float64}

- This function is used within the `exsolve3` function.
- Calculate value from different parameter sets (# of parameter sets: Silicic: 3, Mafic: 2)

# Arguments
- `composition`: Silicic or Mafic
- `s`: represent different set of parammeters, the Silicic case has "a", "b", "c", and the Mafic has "a", "b".
- `x`: Water (H2O)
- `y`: Gas (CO2)
- `z`: Pressure (P)
"""
function dX_dxdydz(
    composition::Union{Silicic,Mafic}, s::String, x::T, y::T, z::T
)::NamedTuple{(:X, :dXdx, :dXdy, :dXdz),NTuple{4,T}} where {T<:Float64}
    @unpack c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 = meltingcurve_dict[composition][s]
    X =
        c1 +
        c2 * x +
        c3 * y +
        c4 * z +
        c5 * x * y +
        c6 * x * z +
        c7 * y * z +
        c8 * x^2 +
        c9 * y^2 +
        c10 * z^2
    dXdx = 100 * (c2 + c5 * y + c6 * z + 2 * c8 * x)
    dXdy = 100 * (c3 + c5 * x + c7 * z + 2 * c9 * y)
    dXdz = 1e-6 * (c4 + c6 * x + c7 * y + 2 * c10 * z)
    return (; X, dXdx, dXdy, dXdz)
end
