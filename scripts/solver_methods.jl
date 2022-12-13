methods = Dict(
    "QNDF" => QNDF(; autodiff=false),
    "FBDF" => FBDF(; autodiff=false),
    "Rodas4" => Rodas4(; autodiff=false),
    "Tsit5" => Tsit5(),
    "KenCarp4" => KenCarp4(; autodiff=false),
    "CVODE_BDF" => CVODE_BDF(),
    "Rosenbrock23" => Rosenbrock23(; autodiff=false),
)
