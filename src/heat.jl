"""
    gas_heat_capacity(X_co2::Float64)

# Arguments
`X_co2`: mole fraction of CO2 in the gas.
"""
function gas_heat_capacity(X_co2::Float64)::Float64
    if X_co2 == 0
        c_g = 0.0
        return c_g
    end

    # Properties of CO2
    m_co2 = 44.01e-3
    c_co2 = 1200.0

    # Properties of H2O
    m_h2o = 18.02e-3
    c_h2o = 3880.0

    # effective molar mass
    m_g = m_h2o*(1-X_co2)+m_co2*X_co2
    c_g = (m_h2o*c_h2o*(1-X_co2)+m_co2*c_co2*X_co2)/m_g

    return c_g
end

"""
    heat_conduction_chamberCH(maxn::Int64, a::Float64, c::Float64, dr::Float64, kappa::Float64, rho::Float64, cp::Float64, Tb::Float64, param_sv)::Float64

# Arguments
`maxn`: number of terms
`a`: radius of magma chamber (m)
`c`: radius of outer shell (m)
`dr`: grid spacing of calculate the heat transfer
`kappa`: thermal diffusivity of host rocks
`rho`: density of the magma
`cp`: specific heat of magma
`Tb`: boundary temperature of the outer shell
"""
function heat_conduction_chamberCH(maxn::Int64, a::Float64, c::Float64, dr::Float64, kappa::Float64, rho::Float64, cp::Float64, Tb::Float64, param_sv::ParamSaved{Float64})::Float64
    storeTime = param_sv.storeTime
    maxTime = param_sv.maxTime
    storeTemp = param_sv.storeTemp
    lengthTime = param_sv.lengthTime
    switch_Tprofile = param_sv.switch_Tprofile
    storeSumk = param_sv.storeSumk
    storeSumk_old = param_sv.storeSumk_old
    storeSumk_2 = param_sv.storeSumk_2
    storeSumk_2_old = param_sv.storeSumk_2_old
    # geometry
    r     = a + dr # radial coordinate (m)

    # time grid
    time = storeTime
    time_index= length(time)

    if maxTime < time[end]
        maxTime = time[end]
        storeSumk_old = copy(storeSumk)
        storeSumk_2_old = copy(storeSumk_2)
    end

    # initial and boundary conditions
    Tr0   = storeTemp
    # pick a time
    current_time    = time[end]

    if time_index == 2
        # first sum over n
        sumn = 0.0
        for n=1:maxn
            # sum over k within first sum over n
            sumk = 0.0
            for k=1:time_index-1
                past_time       = time[k]
                past_time_delta = time[k+1]
                mk    = ((Tr0[k+1] - Tr0[k])/(time[k+1]-time[k]))
                pk    = (Tr0[k] - mk*time[k])
                termk = mk*c^4/(kappa^2*n^4*pi^4)*(exp(-kappa*(n*pi/c)^2*(current_time-past_time))*(1-kappa*(n*pi/c)^2*past_time) + exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta))*(kappa*(n*pi/c)^2*past_time_delta-1)) + pk*c^2/(kappa*n^2*pi^2)*(exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) - exp(-kappa*(n*pi/c)^2*(current_time-past_time)))
                sumk = sumk + termk
            end
            storeSumk[n] = sumk
            termn = n*sin(n*pi*(r-a)/c)*sumk
            sumn  = sumn + termn
        end
        term1 = -4*pi*kappa*a/(2*r*c^2)*(-1)*sumn

        # second sum over n
        sumn = 0.0
        for n=1:maxn
            # sum over k within second sum over n
            sumk = 0.0
            for k=1:time_index-1
                past_time       = time[k]
                past_time_delta = time[k+1]
                termk = exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) - exp(-kappa*(n*pi/c)^2*(current_time-past_time))
                sumk = sumk + termk
            end
            storeSumk_2[n] = sumk
            termn = 1/n*sin(n*pi*(r-a)/c)*cos(n*pi)*sumk
            sumn  = sumn + termn
        end
        term2 = -4*pi*kappa*(a+c)/(2*kappa*pi^2*r)*Tb*(1)*sumn
        lengthTime = 2
        param_sv.lengthTime = lengthTime

    elseif time_index > 2 && time_index > lengthTime
        sumn = 0.0
        for n=1:maxn
            # sum over k within first sum over n
            temp = storeSumk[n]*exp(-kappa*(n*pi/c)^2*(current_time-time[time_index-1]))
            past_time       = time[time_index-1]
            past_time_delta = time[time_index]
            mk    = (Tr0[time_index] - Tr0[time_index-1])/(time[time_index]-time[time_index-1])
            pk    = Tr0[time_index-1] - mk*time[time_index-1]
            termk = mk*c^4/(kappa^2*n^4*pi^4)*(exp(-kappa*(n*pi/c)^2*(current_time-past_time))*(1-kappa*(n*pi/c)^2*past_time) + exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta))*(kappa*(n*pi/c)^2*past_time_delta-1)) + pk*c^2/(kappa*n^2*pi^2)*(exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) - exp(-kappa*(n*pi/c)^2*(current_time-past_time)))
            sumk = temp + termk
            storeSumk_old[n] = storeSumk[n]
            storeSumk[n] = sumk
            termn = n*sin(n*pi*(r-a)/c)*sumk
            sumn  = sumn + termn
        end
        param_sv.storeSumk_old = copy(storeSumk_old)
        term1 = -4*pi*kappa*a/(2*r*c^2)*(-1)*sumn

        # second sum over n
        sumn = 0.0
        for n=1:maxn
            temp = storeSumk_2[n]*exp(-kappa*(n*pi/c)^2*(current_time-time[time_index-1]))
            past_time       = time[time_index-1]
            past_time_delta = time[time_index]
            termk = exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) - exp(-kappa*(n*pi/c)^2*(current_time-past_time))
            sumk = temp + termk
            storeSumk_2_old[n] = storeSumk_2[n]
            storeSumk_2[n] = sumk
            termn = 1/n*sin(n*pi*(r-a)/c)*cos(n*pi)*sumk
            sumn  = sumn + termn
        end
        term2 = -4*pi*kappa*(a+c)/(2*kappa*pi^2*r)*Tb*(1)*sumn
        switch_Tprofile = 0
        param_sv.storeSumk_2_old = copy(storeSumk_2_old)
        param_sv.switch_Tprofile = switch_Tprofile

    elseif time_index > 2 && time_index == lengthTime
        sumn = 0.0
        for n=1:maxn
            # sum over k within first sum over n
            if switch_Tprofile == 0
                temp = storeSumk_old[n]*exp(-kappa*(n*pi/c)^2*(current_time-time[time_index-1]))
                past_time       = time[time_index-1]
                past_time_delta = time[time_index]
                mk    = (Tr0[end] - Tr0[time_index-1])/(time[time_index]-time[time_index-1])
                pk    = Tr0[time_index-1] - mk*time[time_index-1]
                termk = mk*c^4/(kappa^2*n^4*pi^4)*(exp(-kappa*(n*pi/c)^2*(current_time-past_time))*(1-kappa*(n*pi/c)^2*past_time) + exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta))*(kappa*(n*pi/c)^2*past_time_delta-1)) + pk*c^2/(kappa*n^2*pi^2)*(exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) - exp(-kappa*(n*pi/c)^2*(current_time-past_time)))
                sumk = temp + termk
                storeSumk[n] = sumk
            else
                sumk = storeSumk_old[n]
            end
            termn = n*sin(n*pi*(r-a)/c)*sumk
            sumn  = sumn + termn
        end
        term1 = -4*pi*kappa*a/(2*r*c^2)*(-1)*sumn

        # second sum over n
        sumn = 0.0
        for n=1:maxn
            if switch_Tprofile == 0
                temp = storeSumk_2_old[n]*exp(-kappa*(n*pi/c)^2*(current_time-time[end-1]))
                past_time       = time[end-1]
                past_time_delta = time[end]
                termk = exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) - exp(-kappa*(n*pi/c)^2*(current_time-past_time))
                sumk = temp + termk
                storeSumk_2[n] = sumk
            else
                sumk = storeSumk_2_old[n]
            end
            termn = 1/n*sin(n*pi*(r-a)/c)*cos(n*pi)*sumk
            sumn  = sumn + termn
        end
        term2 = -4*pi*kappa*(a+c)/(2*kappa*pi^2*r)*Tb*(1)*sumn
    elseif time_index > 2 && time_index < lengthTime
        sumn = 0.0
        switch_Tprofile = 1
        param_sv.switch_Tprofile = switch_Tprofile
        for n=1:maxn
            termn = n*sin(n*pi*(r-a)/c)*storeSumk_old[n]
            sumn  = sumn + termn
        end
        term1 = -4*pi*kappa*a/(2*r*c^2)*(-1)*sumn
        # second sum over n
        sumn = 0.0
        for n=1:maxn
            termn = 1/n*sin(n*pi*(r-a)/c)*cos(n*pi)*storeSumk_2_old[n]
            sumn  = sumn + termn
        end
        term2 = -4*pi*kappa*(a+c)/(2*kappa*pi^2*r)*Tb*(1)*sumn
    elseif time_index < 2
        term1 = 0.0
        term2 = 0.0
    end

    # third sum over n
    Ta = Tr0[end]
    sumn = 0.0
    for n=1:maxn
        termn = sin(n*pi*(r-a)/c)*exp(-kappa*(n*pi/c)^2*current_time)*((a*(a+c)*Ta - a*(a+c)*Tb)/(n*pi)*(1-cos(n*pi)) + ((a+c)*Tb-a*Ta)/(n*pi)*(-(a+c)*cos(n*pi) + a))
        sumn = sumn + termn
    end

    term3 = 2/(r*c)*sumn
    Trt = term1 + term2 + term3
    lengthTime = time_index
    param_sv.lengthTime = lengthTime
    small_q = -kappa*rho*cp*(Trt-Tr0[end])/dr
    surface_area_chamber = 4*pi*a^2
    Q = small_q*surface_area_chamber
    param_sv.maxTime = maxTime
    return Q
end

"""
    heat_conduction_chamber_profileCH(maxn::Int64, a::Float64, c::Float64, r::Float64, kappa::Float64, Tb::Float64, param_sv)::Float64

# Arguments
`maxn`: number of terms
`a`: radius of magma chamber (m)
`c`: radius of outer shell (m)
`r`: 
`kappa`: thermal diffusivity of host rocks
`Tb`: boundary temperature of the outer shell
"""
function heat_conduction_chamber_profileCH(maxn::Int64, a::Float64, c::Float64, r::Vector{Float64}, kappa::Float64, Tb::Float64, param_sv::ParamSaved{Float64})::Float64
    storeTime = param_sv.storeTime
    storeTemp = param_sv.storeTemp
    lengthTime = param_sv.lengthTime
    switch_Tprofile = param_sv.switch_Tprofile
    storeSumk = param_sv.storeSumk
    storeSumk_old = param_sv.storeSumk_old
    storeSumk_2 = param_sv.storeSumk_2
    storeSumk_2_old = param_sv.storeSumk_2_old

    # temperature history
    time = storeTime
    Tr0  = storeTemp
    time_index= length(time)

    # pick a time
    current_time    = time[end]

    if time_index == 2
        # first sum over n
        sumn = 0.0
        for n=1:maxn
            # sum over k within first sum over n
            sumk = 0.0
            for k=1:time_index-1
                past_time       = time[k]
                past_time_delta = time[k+1]
                mk    = (Tr0[k+1] - Tr0[k])/(time[k+1]-time[k])
                pk    = Tr0[k] - mk*time[k]
                termk = mk*c^4/(kappa^2*n^4*pi^4)*(exp(-kappa*(n*pi/c)^2*(current_time-past_time))*(1-kappa*(n*pi/c)^2*past_time) + exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta))*(kappa*(n*pi/c)^2*past_time_delta-1)) + pk*c^2/(kappa*n^2*pi^2)*(exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) - exp(-kappa*(n*pi/c)^2*(current_time-past_time)))
                sumk = sumk + termk
            end
            termn = n*sin.(n*pi*(r.-a)/c)*sumk
            sumn  = sumn .+ termn
        end
        term1 = -4*pi*kappa*a/(2*r*c^2)*(-1)*sumn

        #% second sum over n
        sumn = 0.0
        for n=1:maxn
            #% sum over k within second sum over n
            sumk = 0.0
            for k=1:time_index-1
                past_time       = time[k]
                past_time_delta = time[k+1]
                termk = exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) - exp(-kappa*(n*pi/c)^2*(current_time-past_time))
                sumk = sumk + termk
            end
            termn = 1/n*sin.(n*pi*(r.-a)/c)*cos(n*pi)*sumk
            sumn  = sumn .+ termn
        end
        term2 = -4*pi*kappa*(a+c)/(2*kappa*pi^2*r)*Tb*(1)*sumn

    elseif time_index > 2 && time_index > lengthTime
        # first sum over n
        sumn = 0.0
        for n=1:maxn
            # sum over k within first sum over n
            temp = storeSumk[n]*exp(-kappa*(n*pi/c)^2*(current_time-time[time_index-1]))
            past_time       = time[time_index-1]
            past_time_delta = time[time_index]
            mk    = (Tr0[time_index] - Tr0[time_index-1])/(time[time_index]-time[time_index-1])
            pk    = Tr0[time_index-1] - mk*time[time_index-1]
            termk = mk*c^4/(kappa^2*n^4*pi^4)*(exp(-kappa*(n*pi/c)^2*(current_time-past_time))*(1-kappa*(n*pi/c)^2*past_time) + exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta))*(kappa*(n*pi/c)^2*past_time_delta-1)) + pk*c^2/(kappa*n^2*pi^2)*(exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) - exp(-kappa*(n*pi/c)^2*(current_time-past_time)))
            sumk = temp + termk
            termn = n*sin.(n*pi*(r.-a)/c)*sumk
            sumn  = sumn .+ termn
        end
        term1 = -4*pi*kappa*a/(2*r*c^2)*(-1)*sumn

        # second sum over n
        sumn = 0.0
        for n=1:maxn
            temp=storeSumk_2[n]*exp(-kappa*(n*pi/c)^2*(current_time-time[time_index-1]))
            past_time       = time[time_index-1]
            past_time_delta = time[time_index]
            termk = exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) - exp(-kappa*(n*pi/c)^2*(current_time-past_time))
            sumk = temp + termk
            termn = 1/n*sin.(n*pi*(r.-a)/c)*cos(n*pi)*sumk
            sumn  = sumn .+ termn
        end
        term2 = -4*pi*kappa*(a+c)/(2*kappa*pi^2*r)*Tb*(1)*sumn
        switch_Tprofile = 0
        param_sv.switch_Tprofile = switch_Tprofile

    elseif time_index > 2 && time_index == lengthTime
        sumn = 0.0
        for n=1:maxn
            # sum over k within first sum over n
            if switch_Tprofile == 0
                temp = storeSumk_old[n]*exp(-kappa*(n*pi/c)^2*(current_time-time[time_index-1]))
                past_time       = time[time_index-1]
                past_time_delta = time[time_index]
                mk    = (Tr0[time_index] - Tr0[time_index-1])/(time[time_index]-time[time_index-1])
                pk    = Tr0[time_index-1] - mk*time[time_index-1]
                termk = mk*c^4/(kappa^2*n^4*pi^4)*(exp(-kappa*(n*pi/c)^2*(current_time-past_time))*(1-kappa*(n*pi/c)^2*past_time) + exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta))*(kappa*(n*pi/c)^2*past_time_delta-1)) + pk*c^2/(kappa*n^2*pi^2)*(exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) - exp(-kappa*(n*pi/c)^2*(current_time-past_time)))
                sumk = temp + termk
            else
                sumk = storeSumk_old[n]
            end
            termn = n*sin.(n*pi*(r.-a)/c)*sumk
            sumn  = sumn .+ termn
        end
        term1 = -4*pi*kappa*a/(2*r*c^2)*(-1)*sumn

        # second sum over n
        sumn = 0.0
        for n=1:maxn
            if switch_Tprofile == 0
                temp=storeSumk_2_old[n]*exp(-kappa*(n*pi/c)^2*(current_time-time[end-1]))
                past_time       = time[end-1]
                past_time_delta = time[end]
                termk = exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) - exp(-kappa*(n*pi/c)^2*(current_time-past_time))
                sumk = temp + termk
            else
                sumk = storeSumk_2_old[n]
            end
            termn = 1/n*sin.(n*pi*(r.-a)/c)*cos(n*pi)*sumk
            sumn  = sumn .+ termn
        end
        term2 = -4*pi*kappa*(a+c)/(2*kappa*pi^2*r)*Tb*(1)*sumn

    elseif time_index > 2 && time_index < lengthTime
        sumn = 0.0
        switch_Tprofile = 1
        param_sv.switch_Tprofile = switch_Tprofile
        sumn = 0.0
        for n=1:maxn
            termn = n*sin.(n*pi*(r.-a)/c)*storeSumk_old[n]
            sumn  = sumn .+ termn
        end
        term1 = -4*pi*kappa*a/(2*r*c^2)*(-1)*sumn

        # second sum over n
        sumn = 0.0
        for n=1:maxn
            termn = 1/n*sin.(n*pi*(r.-a)/c)*cos(n*pi)*storeSumk_2_old[n]
            sumn  = sumn .+ termn
        end
        term2 = -4*pi*kappa*(a+c)/(2*kappa*pi^2*r)*Tb*(1)*sumn

    elseif time_index < 2
        term1 = 0.0
        term2 = 0.0
    end

    # third sum over n
    Ta = Tr0[time_index]
    sumn = 0.0
    for n=1:maxn
        termn =  sin.(n*pi*(r.-a)/c)*exp(-kappa*(n*pi/c)^2*current_time)*((a*(a+c)*Ta - a*(a+c)*Tb)/(n*pi)*(1-cos(n*pi))+((a+c)*Tb -a*Ta)/(n*pi)*(-(a+c)*cos(n*pi) + a))
        sumn = sumn .+ termn
    end
    term3 = 2/(r*c)*sumn
    Trt = term1 .+ term2 .+ term3
    return Trt
end
