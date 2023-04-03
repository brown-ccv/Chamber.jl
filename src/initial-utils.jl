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

meq_silicic(
    Pw::Float64,
    Pc::Float64,
    Temp::Float64,
    b1::Float64,
    b2::Float64,
    b3::Float64,
    b4::Float64,
    b5::Float64,
    b6::Float64,
)::Float64 =
    1e-2 * real(
        (b1 * Complex(Pw)^0.5 + b2 * Pw + b3 * Complex(Pw)^1.5) / Temp +
        b4 * Complex(Pw)^1.5 +
        Pc * (b5 * Complex(Pw)^0.5 + b6 * Pw),
    )
dmeqdT_silicic(Pw::Float64, Temp::Float64, b1::Float64, b2::Float64, b3::Float64)::Float64 =
    1e-2 * real(-1 * (b1 * Complex(Pw)^0.5 + b2 * Pw + b3 * Complex(Pw)^1.5) * Temp^(-2))
dmeqdP_silicic(
    Pw::Float64,
    dPwdP::Float64,
    Pc::Float64,
    dPcdP::Float64,
    Temp::Float64,
    b1::Float64,
    b2::Float64,
    b3::Float64,
    b4::Float64,
    b5::Float64,
    b6::Float64,
)::Float64 =
    1e-8 * real(
        dPwdP * (
            (1 / Temp) * (0.5 * b1 * Complex(Pw)^(-0.5) + b2 + 1.5 * b3 * Complex(Pw)^0.5) +
            1.5 * b4 * Complex(Pw)^0.5 +
            Pc * (0.5 * b5 * Complex(Pw)^(-0.5) + b6)
        ) + dPcdP * (b5 * Complex(Pw)^0.5 + b6 * Pw),
    )
dmeqdXco2_silicic(
    Pw::Float64,
    dPwdXco2::Float64,
    Pc::Float64,
    dPcdXco2::Float64,
    Temp::Float64,
    b1::Float64,
    b2::Float64,
    b3::Float64,
    b4::Float64,
    b5::Float64,
    b6::Float64,
)::Float64 =
    1e-2 * real(
        dPwdXco2 * (
            (1 / Temp) * (0.5 * b1 * Complex(Pw)^(-0.5) + b2 + 1.5 * b3 * Complex(Pw)^0.5) +
            1.5 * b4 * Complex(Pw)^0.5 +
            Pc * (0.5 * b5 * Complex(Pw)^(-0.5) + b6)
        ) + dPcdXco2 * (b5 * Complex(Pw)^0.5 + b6 * Pw),
    )

meq_mafic(
    P::Float64,
    T_C::Float64,
    X_co2::Float64,
    b1::Float64,
    b2::Float64,
    b3::Float64,
    b4::Float64,
    b5::Float64,
    b6::Float64,
    b7::Float64,
    b8::Float64,
    b9::Float64,
    b10::Float64,
)::Float64 =
    1e-2 * real(
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
dmeqdT_mafic(
    P::Float64,
    T_C::Float64,
    X_co2::Float64,
    b2::Float64,
    b5::Float64,
    b6::Float64,
    b8::Float64,
)::Float64 = 1e-2 * real(b2 + b5 * X_co2 + b6 * P + 2 * b8 * T_C)
dmeqdP_mafic(
    P::Float64,
    T_C::Float64,
    X_co2::Float64,
    b4::Float64,
    b6::Float64,
    b7::Float64,
    b10::Float64,
)::Float64 = 1e-8 * real(b4 + b6 * T_C + b7 * X_co2 + 2 * b10 * P)
dmeqdXco2_mafic(
    P::Float64,
    T_C::Float64,
    X_co2::Float64,
    b3::Float64,
    b5::Float64,
    b7::Float64,
    b9::Float64,
)::Float64 = 1e-2 * real(b3 + b5 * T_C + b7 * P + 2 * b9 * X_co2)

C_co2_f(
    Pw::Float64,
    Pc::Float64,
    Temp::Float64,
    c1::Float64,
    c2::Float64,
    c3::Float64,
    c4::Float64,
)::Float64 =
    1e-6 *
    real(Pc * (c1 + c2 * Pw) / Temp + Pc * (c3 * Complex(Pw)^0.5 + c4 * Complex(Pw)^1.5))
dC_co2dT_f(Pw::Float64, Pc::Float64, Temp::Float64, c1::Float64, c2::Float64)::Float64 =
    1e-6 * real(-1 * Pc * (c1 + c2 * Pw) * Complex(Temp)^(-2))
dC_co2dP_f(
    Pw::Float64,
    Pc::Float64,
    Temp::Float64,
    dPwdP::Float64,
    dPcdP::Float64,
    c1::Float64,
    c2::Float64,
    c3::Float64,
    c4::Float64,
)::Float64 =
    1e-12 * real(
        Complex(Temp)^(-1) * (dPcdP * (c1 + c2 * Pw) + Pc * (c2 * dPwdP)) +
        dPcdP * (c3 * Complex(Pw)^0.5 + c4 * Complex(Pw)^1.5) +
        Pc * (0.5 * c3 * Complex(Pw)^(-0.5) * dPwdP + 1.5 * c4 * Complex(Pw)^0.5 * dPwdP),
    )
dC_co2dXco2_f(
    Pw::Float64,
    Pc::Float64,
    Temp::Float64,
    dPwdXco2::Float64,
    dPcdXco2::Float64,
    c1::Float64,
    c2::Float64,
    c3::Float64,
    c4::Float64,
)::Float64 =
    1e-6 * real(
        Complex(Temp)^(-1) * (dPcdXco2 * (c1 + c2 * Pw) + Pc * (c2 * dPwdXco2)) +
        dPcdXco2 * (c3 * Complex(Pw)^0.5 + c4 * Complex(Pw)^1.5) +
        Pc *
        (0.5 * c3 * Complex(Pw)^(-0.5) * dPwdXco2 + 1.5 * c4 * Complex(Pw)^0.5 * dPwdXco2),
    )

# For function exsolve
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

# For function exsolve
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
function water(
    composition::Silicic, p::Float64, t::Float64, x::Float64, c::Float64
)::Float64
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

function water(composition::Mafic, p::Float64, t::Float64, x::Float64, c::Float64)::Float64
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
function dwater_dx(composition::Silicic, p::Float64, t::Float64, x::Float64)::Float64
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

function dwater_dx(composition::Mafic, p::Float64, t::Float64, x::Float64)::Float64
    @unpack h3, h5, h7, h9 = Exsolve3Mafic()
    return h3 + h5 * t + h7 * p + 2 * h9 * x
end

# For function exsolve3, finding X_CO2
function solve_NR(
    f, f_prime, errorTol::Float64, count_max::Float64, Xc_initial::Float64
)::Float64
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
a_f(x::Float64, y::Float64, z::Float64)::Float64 = 0.36 - 0.02x - 0.06y + 8.6e-4z + 0.0024x*y + 6.27e-5x*z + 3.57e-5y*z - 0.0026x^2 + 0.003y^2 - 1.16e-6z^2
dadx_f(x::Float64, y::Float64, z::Float64)::Float64 = 100(-0.02 + 0.0024y + 6.27e-5z - 2*0.0026x)
dady_f(x::Float64, y::Float64, z::Float64)::Float64 = 100(-0.06 + 0.0024x + 3.57e-5z + 2*0.003y)
dadz_f(x::Float64, y::Float64, z::Float64)::Float64 = 1e-6(8.6e-4 + 6.27e-5x + 3.57e-5y - 2*1.16e-6z)

b_f(x::Float64, y::Float64, z::Float64)::Float64 = 0.0071 + 0.0049x + 0.0043y - 4.08e-5z - 7.85e-4x*y - 1.3e-5x*z + 3.97e-6y*z + 6.29e-4x^2 - 0.0025y^2 + 8.51e-8z^2
dbdx_f(x::Float64, y::Float64, z::Float64)::Float64 = 100(0.0049 - 7.85e-4y - 1.3e-5z + 2*6.29e-4x)
dbdy_f(x::Float64, y::Float64, z::Float64)::Float64 = 100(0.0043 - 7.85e-4x + 3.97e-6z - 2*0.0025y)
dbdz_f(x::Float64, y::Float64, z::Float64)::Float64 = 1e-6(-4.08e-5 - 1.3e-5x + 3.97e-6y + 2*8.51e-8z)

c_f(x::Float64, y::Float64, z::Float64)::Float64 = 863.09 - 36.9x + 48.81y - 0.17z - 1.52x*y - 0.04x*z - 0.04y*z + 4.57x^2 - 7.79y^2 + 4.65e-4z^2
dcdx_f(x::Float64, y::Float64, z::Float64)::Float64 = 100(-36.9 - 1.52y - 0.04z + 2*4.57x)
dcdy_f(x::Float64, y::Float64, z::Float64)::Float64 = 100(48.81 - 1.52x - 0.04z - 2*7.79y)
dcdz_f(x::Float64, y::Float64, z::Float64)::Float64 = 1e-6(-0.17 - 0.04x - 0.04y + 2*4.65e-4z)

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
