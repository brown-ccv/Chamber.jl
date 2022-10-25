# For function eos_g
rho_g(P, T)     = real(1e3*(-112.528*Complex(T-273.15)^-0.381 + 127.811*Complex(P*1e-5)^-1.135 + 112.04*Complex(T-273.15)^-0.411*Complex(P*1e-5)^0.033))
drho_g_dP(P, T) = real(1e-2*((-1.135)*127.811*Complex(P*1e-5)^-2.135 + 0.033*112.04*Complex(T-273.15)^-0.411*Complex(P*1e-5)^-0.967))
drho_g_dT(P, T) = real(1e3*((-0.381)*(-112.528)*Complex(T-273.15)^-1.381 + (-0.411)*112.04*Complex(T-273.15)^-1.411*Complex(P*1e-5)^0.033))

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