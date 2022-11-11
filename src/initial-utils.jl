# For function exsolve
struct ExsolveSilicic{T}
    b1::T
    b2::T
    b3::T
    b4::T
    b5::T
    b6::T
    ExsolveSilicic() = new{Float64}(354.94,
                                    9.623,
                                    -1.5223,
                                    1.2439e-3,
                                    -1.084e-4,
                                    -1.362e-5)
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
    ExsolveMafic() = new{Float64}(2.99622526644026,
                                  0.00322422830627781,
                                  -9.1389095360385,
                                  0.0336065247530767,
                                  0.00747236662935722,
                                  -0.0000150329805347769,
                                  -0.01233608521548,
                                  -4.14842647942619e-6,
                                  -0.655454303068124,
                                  -7.35270395041104e-6)
end

struct Co2PartitionCoeff{T}
    c1::T
    c2::T
    c3::T
    c4::T
    Co2PartitionCoeff() = new{Float64}(5668,
                                       -55.99,
                                       0.4133,
                                       2.041e-3)
end

meq_silicic(Pw, Pc, Temp, b1, b2, b3, b4, b5, b6)                           = 1e-2*real((b1*Complex(Pw)^0.5+b2*Pw+b3*Complex(Pw)^1.5)/Temp+b4*Complex(Pw)^1.5+Pc*(b5*Complex(Pw)^0.5+b6*Pw))
dmeqdT_silicic(Pw, Temp, b1, b2, b3)                                        = 1e-2*real(-1*(b1*Complex(Pw)^0.5+b2*Pw+b3*Complex(Pw)^1.5)*Temp^(-2))
dmeqdP_silicic(Pw, dPwdP, Pc, dPcdP, Temp, b1, b2, b3, b4, b5, b6)          = 1e-8*real(dPwdP*((1/Temp)*(0.5*b1*Complex(Pw)^(-0.5)+b2+1.5*b3*Complex(Pw)^0.5)+1.5*b4*Complex(Pw)^0.5+Pc*(0.5*b5*Complex(Pw)^(-0.5)+b6))+dPcdP*(b5*Complex(Pw)^0.5+b6*Pw))
dmeqdXco2_silicic(Pw, dPwdXco2, Pc, dPcdXco2, Temp, b1, b2, b3, b4, b5, b6) = 1e-2*real(dPwdXco2*((1/Temp)*(0.5*b1*Complex(Pw)^(-0.5)+b2+1.5*b3*Complex(Pw)^0.5)+1.5*b4*Complex(Pw)^0.5+Pc*(0.5*b5*Complex(Pw)^(-0.5)+b6))+dPcdXco2*(b5*Complex(Pw)^0.5+b6*Pw))

meq_mafic(P, T_C, X_co2, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10) = 1e-2*real(b1+b2*T_C+b3*X_co2+b4*P+b5*T_C*X_co2+b6*T_C*P+b7*X_co2*P+b8*T_C^2+b9*Complex(X_co2)^2+b10*Complex(P)^2)
dmeqdT_mafic(P, T_C, X_co2, b2, b5, b6, b8)                       = 1e-2*real(b2+b5*X_co2+b6*P+2*b8*T_C)
dmeqdP_mafic(P, T_C, X_co2, b4, b6, b7, b10)                      = 1e-8*real(b4+b6*T_C+b7*X_co2+2*b10*P)
dmeqdXco2_mafic(P, T_C, X_co2, b3, b5, b7, b9)                    = 1e-2*real(b3+b5*T_C+b7*P+2*b9*X_co2)

C_co2_f(Pw, Pc, Temp, c1, c2, c3, c4)                             = 1e-6*real(Pc*(c1+c2*Pw)/Temp+Pc*(c3*Complex(Pw)^0.5+c4*Complex(Pw)^1.5))
dC_co2dT_f(Pw, Pc, Temp, c1, c2)                                  = 1e-6*real(-1*Pc*(c1+c2*Pw)*Complex(Temp)^(-2))
dC_co2dP_f(Pw, Pc, Temp, dPwdP, dPcdP, c1, c2, c3, c4)            = 1e-12*real(Complex(Temp)^(-1)*(dPcdP*(c1+c2*Pw)+Pc*(c2*dPwdP))+dPcdP*(c3*Complex(Pw)^0.5+c4*Complex(Pw)^1.5)+Pc*(0.5*c3*Complex(Pw)^(-0.5)*dPwdP+1.5*c4*Complex(Pw)^0.5*dPwdP))
dC_co2dXco2_f(Pw, Pc, Temp, dPwdXco2, dPcdXco2, c1, c2, c3, c4)   = 1e-6*real(Complex(Temp)^(-1)*(dPcdXco2*(c1+c2*Pw)+Pc*(c2*dPwdXco2))+dPcdXco2*(c3*Complex(Pw)^0.5+c4*Complex(Pw)^1.5)+Pc*(0.5*c3*Complex(Pw)^(-0.5)*dPwdXco2+1.5*c4*Complex(Pw)^0.5*dPwdXco2))

function build_meq_silicic(Pw::T, Pc::T, Temp::T, dPwdP::T, dPcdP::T, dPwdXco2::T, dPcdXco2::T)::Vector{T} where T<:Float64
    @unpack b1, b2, b3, b4, b5, b6 = ExsolveSilicic()
    meq       = meq_silicic(Pw, Pc, Temp, b1, b2, b3, b4, b5, b6)
    dmeqdT    = dmeqdT_silicic(Pw, Temp, b1, b2, b3)
    dmeqdP    = dmeqdP_silicic(Pw, dPwdP, Pc, dPcdP, Temp, b1, b2, b3, b4, b5, b6)
    dmeqdXco2 = dmeqdXco2_silicic(Pw, dPwdXco2, Pc, dPcdXco2, Temp, b1, b2, b3, b4, b5, b6)
    return [meq, dmeqdT, dmeqdP, dmeqdXco2]
end

function build_meq_mafic(P::T, Temp::T, X_co2::T)::Vector{T} where T<:Float64
    T_C = Temp-273.15
    @unpack b1, b2, b3, b4, b5, b6, b7, b8, b9, b10 = ExsolveMafic()
    meq       = meq_mafic(P, T_C, X_co2, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
    dmeqdT    = dmeqdT_mafic(P, T_C, X_co2, b2, b5, b6, b8)
    dmeqdP    = dmeqdP_mafic(P, T_C, X_co2, b4, b6, b7, b10)
    dmeqdXco2 = dmeqdXco2_mafic(P, T_C, X_co2, b3, b5, b7, b9)
    return [meq, dmeqdT, dmeqdP, dmeqdXco2]
end

function build_co2(Pw::T, Pc::T, Temp::T, dPwdP::T, dPcdP::T, dPwdXco2::T, dPcdXco2::T)::Vector{T} where T<:Float64
    @unpack c1, c2, c3, c4 = Co2PartitionCoeff()
    C_co2       = C_co2_f(Pw, Pc, Temp, c1, c2, c3, c4)
    dC_co2dT    = dC_co2dT_f(Pw, Pc, Temp, c1, c2)
    dC_co2dP    = dC_co2dP_f(Pw, Pc, Temp, dPwdP, dPcdP, c1, c2, c3, c4)
    dC_co2dXco2 = dC_co2dXco2_f(Pw, Pc, Temp, dPwdXco2, dPcdXco2, c1, c2, c3, c4)
    return [C_co2, dC_co2dT, dC_co2dP, dC_co2dXco2]
end

# For function eos_g
rho_g(P::Float64, T::Float64)::Float64     = real(1e3*(-112.528*Complex(T-273.15)^-0.381 + 127.811*Complex(P*1e-5)^-1.135 + 112.04*Complex(T-273.15)^-0.411*Complex(P*1e-5)^0.033))
drho_g_dP(P::Float64, T::Float64)::Float64 = real(1e-2*((-1.135)*127.811*Complex(P*1e-5)^-2.135 + 0.033*112.04*Complex(T-273.15)^-0.411*Complex(P*1e-5)^-0.967))
drho_g_dT(P::Float64, T::Float64)::Float64 = real(1e3*((-0.381)*(-112.528)*Complex(T-273.15)^-1.381 + (-0.411)*112.04*Complex(T-273.15)^-1.411*Complex(P*1e-5)^0.033))

struct EosG{T}
    rho_g::T
    drho_g_dP::T
    drho_g_dT::T
end
EosG(P, T) = EosG{typeof(P)}(rho_g(P, T),
drho_g_dP(P, T),
drho_g_dT(P, T))

struct EosG_RhoG{T}
    rho_g::T
end
EosG_RhoG(P, T) = EosG_RhoG{typeof(P)}(rho_g(P, T))

