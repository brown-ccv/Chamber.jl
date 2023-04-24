using Parameters

struct Silicic end
struct Mafic end

function get_timestamp()::String
    return Dates.format(now(), "YYmmddHHMM")
end

## Settings for ODE Solver Calculation
@with_kw struct OdeSetting{T}
    reltol::T = 1e-8    # Relative tolerance
    abstol::T = 1e-8    # Absolute tolerance
    first_step::T = 1e5 # Initial step size
    max_step::T = 1e7   # Maximum allowed step size
end

@with_kw struct ConstantValues{T}
    T_surface::T = 273.0     # surface temperature (K)
    T_gradient::T = 32 / 1e3 # thermal gradient (K/m)
    grav_acc::T = 9.81       # gravitational acceleration (m/s2)
    mm_co2::T = 44.01e-3     # molecular mass of CO2
    mm_h2o::T = 18.02e-3     # molecular mass of H2O
end

## Rheology of the crust
struct RheolComposition{T}
    rho_m0::T # initial melt density (kg/m^3)
    rho_x0::T # initial crystal density (kg/m^3)
    c_m::T    # specific heat of melt (J kg^-1 K-1)
    c_x::T    # specific heat of crystals (J kg^-1 K-1)
    L_m::T    # latent heat of melting (J kg^-1) value used by Caricchi and Blundy (2015)
end
silicic = RheolComposition{Float64}(2400, 2600, 1200, 1200, 290e3)
mafic = RheolComposition{Float64}(2420, 2900, 1142, 1160, 470e3)
rheol_composition_dict = Dict("Silicic" => silicic, "Mafic" => mafic)

struct RheolNew{T}
    A::T  # material-dependent constant for viscosity law (Pa s)
    B::T  # molar gas constant (J/mol/K)
    G::T  # activation energy for creep (J/mol)
end
new = RheolNew{Float64}(4.25e7, 8.31, 141e3)
struct RheolOld{T}
    nn::T  # powerlaw exponent
    G::T   # material dependent constant (Pa^-n s^-1)
    M::T   # activation energy for creep (J/mol)
    AA::T  # gas constant = 2e-4*(1e6)^-(nn)
end
old = RheolOld{Float64}((1.9), 141e3, 8.31, 2e-4 * (1e6)^-(1.9))
rheol_dict = Dict("new" => new, "old" => old)

## Parameters
@with_kw mutable struct Param{T}
    composition::Union{Silicic,Mafic}
    rheol::String
    fluxing::Bool = false
    single_eruption::Bool = false
    beta_m::T = 1e10
    beta_x::T = 1e10
    alpha_m::T = 1e-5
    alpha_x::T = 1e-5
    L_e::T = 610e3
    mm_co2::T = 44.01e-3
    mm_h2o::T = 18.02e-3
    alpha_r::T = 1e-5
    beta_r::T = 1e10
    kappa::T = 1e-6
    rho_r::T = 2750
    c_r::T = 1200
    ini_eps_x::T = 0.15
    maxn::Int64 = 10000
    GLQ_n::Int64 = 64
    Q_out_old::T = 0
    DP_crit::T = 20e6
    P_lit::T = 0
    P_lit_0::T = 0
    dP_lit_dt::T = 0
    dP_lit_dt_0::T = 0
    P_lit_drop_max::T = 9e6
    Mdot_in_pass::T = 0
    Mdot_out_pass::T = 10000
    XCO2_in::T = 0.8
    rho_m0::T
    rho_x0::T
    c_m::T
    c_x::T
    L_m::T
    A::T = 0
    B::T = 0
    nn::T = 0
    AA::T = 0
    G::T = 0
    M::T = 0
    Tb::T = 0
    tot_h2o_frac_in::T = 0
    tot_co2_frac_in::T = 0
    T_in::T = 0
end

## eruption/cooling_module/viscous_relaxation control
@with_kw mutable struct SW{T}
    heat_cond::T = 1
    visc_relax::T = 1
    eruption::T = 0
end

@with_kw mutable struct ParamICFinder{T}
    max_count::Int64 = 100
    min_eps_g::T = 1e-10
    eps_g_guess_ini::T = 1e-2
    X_co2_guess_ini::T = 0.2
    fraction::T = 0.2
    delta_X_co2::T = 1e-2
    Tol::T = 0.0
end

@with_kw mutable struct ParamSaved{T}
    maxTime::Number = 0
    lengthTime::Int64 = 0
    switch_Tprofile::Int8 = 0
    phase::T = 0
    storeTime::Vector{T} = []
    storeTemp::Vector{T} = []
    storeSumk::Vector{T} = []
    storeSumk_2::Vector{T} = []
    storeSumk_old::Vector{T} = []
    storeSumk_2_old::Vector{T} = []
end

@with_kw mutable struct EruptSaved{T}
    time::Vector{T} = []
    rho::T = 0
    duration::Vector{T} = []
    mass::Vector{T} = []
    volume::Vector{T} = []
end

function rho_f(;
    eps_m::T, eps_g::T, eps_x::T, rho_m::T, rho_g::T, rho_x::T
)::T where {T<:Float64}
    return eps_m * rho_m + eps_g * rho_g + eps_x * rho_x
end
function drho_dX_f(;
    eps_m::T,
    eps_g::T,
    eps_x::T,
    drho_m_dX::T,
    drho_g_dX::T,
    drho_x_dX::T,
    rho_x::T,
    rho_m::T,
    deps_x_dX::T,
)::T where {T<:Float64}
    return eps_m * drho_m_dX +
           eps_g * drho_g_dX +
           eps_x * drho_x_dX +
           (rho_x - rho_m) * deps_x_dX
end

function rc_f(;
    rho_x::T, eps_x::T, c_x::T, rho_m::T, eps_m::T, c_m::T, rho_g::T, eps_g::T, c_g::T
)::T where {T<:Float64}
    return rho_x * eps_x * c_x + rho_m * eps_m * c_m + rho_g * eps_g * c_g
end

function drc_dX_f(;
    eps_x::T,
    c_x::T,
    drho_x_dX::T,
    eps_g::T,
    c_g::T,
    drho_g_dX::T,
    eps_m::T,
    c_m::T,
    drho_m_dX::T,
    rho_x::T,
    rho_m::T,
    deps_x_dX::T,
)::T where {T<:Float64}
    return eps_x * c_x * drho_x_dX +
           eps_g * c_g * drho_g_dX +
           eps_m * c_m * drho_m_dX +
           (rho_x * c_x - rho_m * c_m) * deps_x_dX
end

function build_rho_rc(
    eps_m::T,
    eps_g::T,
    eps_x::T,
    rho_m::T,
    rho_g::T,
    rho_x::T,
    drho_m_dP::T,
    drho_g_dP::T,
    drho_x_dP::T,
    drho_m_dT::T,
    drho_g_dT::T,
    drho_x_dT::T,
    c_x::T,
    c_m::T,
    c_g::T,
    deps_x_dP::T,
    deps_x_dT::T,
)::Vector{T} where {T<:Float64}
    rho = rho_f(;
        eps_m=eps_m, eps_g=eps_g, eps_x=eps_x, rho_m=rho_m, rho_g=rho_g, rho_x=rho_x
    )
    drho_dP = drho_dX_f(;
        eps_m=eps_m,
        eps_g=eps_g,
        eps_x=eps_x,
        drho_m_dX=drho_m_dP,
        drho_g_dX=drho_g_dP,
        drho_x_dX=drho_x_dP,
        rho_x=rho_x,
        rho_m=rho_m,
        deps_x_dX=deps_x_dP,
    )
    drho_dT = drho_dX_f(;
        eps_m=eps_m,
        eps_g=eps_g,
        eps_x=eps_x,
        drho_m_dX=drho_m_dT,
        drho_g_dX=drho_g_dT,
        drho_x_dX=drho_x_dT,
        rho_x=rho_x,
        rho_m=rho_m,
        deps_x_dX=deps_x_dT,
    )
    drho_deps_g = -rho_m + rho_g

    # computing the product of density and specific heat for the mixture and
    # its derivatives
    rc = rc_f(;
        rho_x=rho_x,
        eps_x=eps_x,
        c_x=c_x,
        rho_m=rho_m,
        eps_m=eps_m,
        c_m=c_m,
        rho_g=rho_g,
        eps_g=eps_g,
        c_g=c_g,
    )
    drc_dP = drc_dX_f(;
        eps_x=eps_x,
        c_x=c_x,
        drho_x_dX=drho_x_dP,
        eps_g=eps_g,
        c_g=c_g,
        drho_g_dX=drho_g_dP,
        eps_m=eps_m,
        c_m=c_m,
        drho_m_dX=drho_m_dP,
        rho_x=rho_x,
        rho_m=rho_m,
        deps_x_dX=deps_x_dP,
    )
    drc_dT = drc_dX_f(;
        eps_x=eps_x,
        c_x=c_x,
        drho_x_dX=drho_x_dT,
        eps_g=eps_g,
        c_g=c_g,
        drho_g_dX=drho_g_dT,
        eps_m=eps_m,
        c_m=c_m,
        drho_m_dX=drho_m_dT,
        rho_x=rho_x,
        rho_m=rho_m,
        deps_x_dX=deps_x_dT,
    )
    return [rho, drho_dP, drho_dT, drho_deps_g, rc, drc_dP, drc_dT]
end

rho_0_f(
    eps_g0::Float64, eps_x0::Float64, rho_g0::Float64, rho_m0::Float64, rho_x0::Float64
)::Float64 = (1 - eps_g0 - eps_x0) * rho_m0 + eps_g0 * rho_g0 + eps_x0 * rho_x0

function build_mdot_in(
    fluxing::Bool, rho_m0::Float64, log_vfr::Float64, P_0::Float64, T_in::Float64
)::Float64
    if ~fluxing
        range_vfr = 10^log_vfr   # volume flow rate (km3/yr)  
        mdot_in = rho_m0 * range_vfr * 1e9 / (3600 * 24 * 365)
    else
        log_vfr = -4.3         # log volume flow rate (km3/yr)
        range_vfr = 10^log_vfr   # volume flow rate (km3/yr)
        rho_g_in = eos_g_rho_g(P_0, T_in)
        mdot_in = rho_g_in * range_vfr * 1e9 / (3600 * 24 * 365)
    end
    return mdot_in
end

function compute_dXdP_dXdT(
    u::Float64, param::Param, var::String
)::Tuple{Float64,Float64,Float64}
    α, β = "alpha_$var", "beta_$var"
    return (u, u / getproperty(param, Symbol(β)), -u * getproperty(param, Symbol(α)))
end

"""
    record_erupt_start(time::T, eps_g::T, eps_x::T, rho_m::T, rho_x::T, rho_g::T, erupt_saved::EruptSaved{T}) where {T<:Float64}

Record time, density of eruptions to `EruptSaved`
- This function is used within the `affect!` function.
"""
function record_erupt_start(
    time::T, eps_g::T, eps_x::T, rho_m::T, rho_x::T, rho_g::T, erupt_saved::EruptSaved{T}
)::Nothing where {T<:Float64}
    push!(erupt_saved.time, time)
    erupt_saved.rho = (1 - eps_g - eps_x) * rho_m + eps_g * rho_g + eps_x * rho_x
    return nothing
end

"""
    record_erupt_end(time::T, erupt_saved::EruptSaved{T}, param::Param{T}) where {T<:Float64}

Record duration, mass and volume of eruptions to `EruptSaved`
- This function is used within the `affect!` function.
"""
function record_erupt_end(
    time::T, erupt_saved::EruptSaved{T}, param::Param{T}
)::Nothing where {T<:Float64}
    duration = time - erupt_saved.time[end]
    mass = duration * param.Mdot_out_pass
    volume = mass / erupt_saved.rho
    push!(erupt_saved.duration, duration)
    push!(erupt_saved.mass, mass)
    push!(erupt_saved.volume, volume)
    return nothing
end
