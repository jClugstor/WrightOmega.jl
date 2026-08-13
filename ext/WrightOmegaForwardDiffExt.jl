module WrightOmegaForwardDiffExt

# import (not using): this file adds a method to wrightomega. Without it, Dual (<: Real)
# falls into the generic wrightomega(::Real) method, which fails at Float64(::Dual).
import WrightOmega: wrightomega
using ForwardDiff: Dual, value, partials

# d/dx wrightomega(x) = w/(1 + w) with w = wrightomega(x)  (TOMS 917, eq. 3).
# value(d) may itself be a Dual (higher-order derivatives); the recursion handles it.
function wrightomega(d::Dual{T}) where {T}
    w = wrightomega(value(d))
    return Dual{T}(w, partials(d) * (w / (1 + w)))
end

end
