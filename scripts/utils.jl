function get_timestamp()
    return Dates.format(now(), "YYmmddHHMM")
end

## Settings for ODE Solver Calculation
struct OdeSetting{T}
    reltol::T     # Relative tolerance
    abstol::T     # Absolute tolerance
    first_step::T # Initial step size
    max_step::T   # Maximum allowed step size
    end

function makeOdeSetting(;reltol=1e-8, abstol=1e-8, first_step=1e5, max_step=1e7)
    return OdeSetting(reltol, abstol, first_step, max_step)
    end

OdeSetting() = OdeSetting{Float64}(1e-8, 1e-8, 1e5, 1e7)


## Rheology of the crust
struct RheolComposition{T}
    rho_m0::T # initial melt density (kg/m^3)
    rho_x0::T # initial crystal density (kg/m^3)
    c_m::T    # specific heat of melt (J kg^-1 K-1)
    c_x::T    # specific heat of crystals (J kg^-1 K-1)
    L_m::T    # latent heat of melting (J kg^-1) value used by Caricchi and Blundy (2015)
end
silicic = RheolComposition{Float64}(2400,
                                    2600,
                                    1200,
                                    1200,
                                    290e3)

mafic =   RheolComposition{Float64}(2420,
                                    2900,
                                    1142,
                                    1160,
                                    470e3)
rheol_composition_dict = Dict("silicic" => silicic, "mafic"  => mafic)

struct RheolNew{T}
    A::T  # material-dependent constant for viscosity law (Pa s)
    B::T  # molar gas constant (J/mol/K)
    G::T  # activation energy for creep (J/mol)
end
new = RheolNew{Float64}(
    4.25e7,
    8.31,
    141e3
)
struct RheolOld{T}
    nn::T  # powerlaw exponent
    G::T   # material dependent constant (Pa^-n s^-1)
    M::T   # activation energy for creep (J/mol)
    AA::T  # gas constant = 2e-4*(1e6)^-(nn)
end
old = RheolOld{Float64}(
    (1.9),
    141e3,
    8.31,
    2e-4*(1e6)^-(1.9)
)
rheol_dict = Dict("new" => new, "old" => old)

## eruption/cooling_module/viscous_relaxation control
mutable struct SW{T}
    heat_cond::T
    visc_relax::T
    eruption::T
end
s = SW{Int8}(1,
             1,
             0)

# Parameters
struct Param{T}
    fluxing::Bool
    single_eruption::Bool
    beta_m::T
    beta_x::T
    alpha_m::T
    alpha_x::T
    L_e::T
    mm_co2::T
    mm_h2o::T
    alpha_r::T
    beta_r::T
    kappa::T
    rho_r::T
    c_r::T
    maxn::Int64
    GLQ_n::Int64
    Q_out_old::T
    dP_lit_dt::T
    dP_lit_dt_0::T
    P_lit_drop_max::T
    Mdot_out_pass::T
    end

function Param{T}(fluxing::Bool=true,
    single_eruption::Bool=true,
    beta_m::T=0,
    beta_x::T=0,
    alpha_m::T=0,
    alpha_x::T=0,
    L_e::T=0,
    mm_co2::T=0,
    mm_h2o::T=0,
    alpha_r::T=0,
    beta_r::T =0,
    kappa::T  =0,
    rho_r::T  =0,
    c_r::T=0,
    maxn::Int64=0,
    GLQ_n::Int64=0,
    Q_out_old::T=0,
    dP_lit_dt::T=0,
    dP_lit_dt_0::T=0,
    P_lit_drop_max::T=0,
    Mdot_out_pass::T=0) where T
    return Param{T}(rheol,
    fluxing,
    single_eruption,
    beta_m,
    beta_x,
    alpha_m,
    alpha_x,
    L_e,
    mm_co2,
    mm_h2o,
    alpha_r,
    beta_r,
    kappa,
    rho_r,
    c_r,
    maxn,
    GLQ_n,
    Q_out_old,
    dP_lit_dt,
    dP_lit_dt_0,
    P_lit_drop_max,
    Mdot_out_pass)
end
p =  Param{Float64}(false,
                    false,
                    1e10,
                    1e10,
                    1e-5,
                    1e-5,
                    610e3,
                    44.01e-3,
                    18.02e-3,
                    1e-5,
                    1e10,
                    1e-6,
                    2750,
                    1200,
                    10000,
                    64,
                    0,
                    0,
                    0,
                    9e6,
                    10000)

mutable struct ParamSaved
    maxTime::Number
    lengthTime::Int64
    switch_Tprofile::Int8
end
ps = ParamSaved(0,
                0,
                0)

mutable struct ParamICFinder
    max_count::Int64
    min_eps_g::Float64
    eps_g_guess_ini::Float64
    X_co2_guess_ini::Float64
    fraction::Float64
    delta_X_co2::Float64
end

pic = ParamICFinder(100,
                    1e-10,
                    1e-2,
                    0.2,
                    0.2,
                    1e-2)