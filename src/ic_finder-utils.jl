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
