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

