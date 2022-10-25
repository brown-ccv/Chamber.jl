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