"""
    gas_heat_capacity(X_co2::Float64)


X_co2
"""
function gas_heat_capacity(X_co2::Float64)
    # %Properties of CO2
    m_co2 = 44.01e-3
    c_co2 = 1200

    # %Properties of H2O
    m_h2o = 18.02e-3
    c_h2o = 3880

    # %effective molar mass
    m_g = m_h2o*(1-X_co2)+m_co2*X_co2
    c_g         = (m_h2o*c_h2o*(1-X_co2)+m_co2*c_co2*X_co2)/m_g
    dc_g_dX_co2 = ((-m_h2o*c_h2o+m_co2*c_co2)+ c_g*(m_h2o-m_co2))/m_g

    if X_co2 == 0
        c_g=0
        dc_g_dX_co2=0
    end
    return [c_g, dc_g_dX_co2]
end