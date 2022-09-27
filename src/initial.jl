"""
    eos_g(P::Number,T::Number)

Initializes some constants. Returns some variables    

# Arguments
-`P`: represents pressure
-`T`: represents the temperature in some units
"""
function eos_g(P::Number,T::Number)
    # parametrization of redlich kwong taken from Huber et al. 2010
    rho_g     = -112.528*Complex(T-273.15)^-0.381 + 127.811*Complex(P*1e-5)^-1.135 + 112.04*Complex(T-273.15)^-0.411*Complex(P*1e-5)^0.033
    drho_g_dP = (-1.135)*127.811*Complex(P*1e-5)^-2.135 + 0.033*112.04*Complex(T-273.15)^-0.411*Complex(P*1e-5)^-0.967
    drho_g_dT = (-0.381)*(-112.528)*Complex(T-273.15)^-1.381 + (-0.411)*112.04*Complex(T-273.15)^-1.411*Complex(P*1e-5)^0.033
    
    rho_g     = real(rho_g*1e3)
    drho_g_dP = real(drho_g_dP*1e-2)
    drho_g_dT = real(drho_g_dT*1e3)
    return Dict(["rho_g"=>rho_g, "drho_g_dP"=>drho_g_dP, "drho_g_dT"=>drho_g_dT])
end