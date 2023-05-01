# find_liq silicic/mafic
h2o, co2, P, ini_eps_x = 0.04, 0.001, 2.1582e8, 0.15
Tl_silicic = 1046.7118332590333
Tl_mafic = 1395.37583723001

# crystal_fraction silicic
T_s, P_s, mH2O_s, mCO2_s = 1046.7118332590333, 2.1582e8, 0.04, 0.001
eps_x_s, deps_x_dP_s, deps_x_dT_s, deps_x_deps_g_s, deps_x_dmco2_t_s, deps_x_dmh2o_t_s = 0.1490308529287438,
-1.0286490935201157e-10, -0.006447776797283597, -1.0, 18.01856777976547,
-13.662995655454651

# crystal_fraction mafic
T_m, P_m, mH2O_m, mCO2_m = 1395.37583723001, 2.1582e8, 0.04, 0.001
eps_x_m, deps_x_dP_m, deps_x_dT_m, deps_x_deps_g_m, deps_x_dmco2_t_m, deps_x_dmh2o_t_m = 0.14913102910967346,
-8.339551728771927e-11, -0.005793139268849489, -1.0, 3.5427927565201927,
16.734989206844432

# parameter_melting_curve silicic
a_s, dadx_s, dady_s, dadz_s, b_s, dbdx_s, dbdy_s, dbdz_s, c_s, dcdx_s, dcdy_s, dcdz_s = (
    0.4912593243499999,
    -0.6673686,
    -5.2193226,
    3.618413000000001e-10,
    0.00234825788664,
    0.21438750000000004,
    0.51204054,
    -4.583466000000008e-12,
    846.285605276,
    -4516.871999999999,
    4010.082,
    2.9072600000000003e-8,
)

# parameter_melting_curve mafic
a_m, dadx_m, dady_m, dadz_m, b_m, dbdx_m, dbdy_m, dbdz_m = (
    -0.009062540132357898,
    0.5874371459003171,
    -0.028483003606881403,
    5.639379444873681e-12,
    10.617411308024828,
    -691.1139598893698,
    31.171886257659274,
    -5.653301393543611e-9,
)
