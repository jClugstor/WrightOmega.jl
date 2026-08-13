module WrightOmega

export wrightomega

include("real.jl")     # Fukushima (2020) piecewise minimax rationals; the fast path
include("complex.jl")  # TOMS Algorithm 917 over the complex plane

end
