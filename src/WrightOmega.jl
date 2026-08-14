module WrightOmega

export wrightomega

include("real.jl")     # Fukushima (2020) piecewise minimax rationals; the fast path
include("complex.jl")  # TOMS Algorithm 917 over the complex plane

"""
    ω

Unicode alias (`\\omega<tab>`) for [`wrightomega`](@ref). Deliberately not exported --
`ω` is a common local variable name -- so use it as `WrightOmega.ω` or bring it in
explicitly with `using WrightOmega: ω`.
"""
const ω = wrightomega

# `public` is a 1.11 keyword and does not parse on 1.10, hence the Expr form.
VERSION >= v"1.11" && eval(Expr(:public, :ω))

end
