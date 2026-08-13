module WrightOmegaMooncakeExt

using WrightOmega: wrightomega
using Mooncake: Mooncake, CoDual, Dual, NoRData, MinimalCtx, @is_primitive,
                zero_fcodual, primal, tangent

# MinimalCtx: a hand-written rule is required for correctness, not just speed --
# the real kernel is @fastmath and Mooncake should never trace into it.
# With no mode argument this registers the signature as primitive in both
# forward and reverse mode.
@is_primitive MinimalCtx Tuple{typeof(wrightomega),Float64}

# d/dx wrightomega(x) = w/(1 + w) with w = wrightomega(x)  (TOMS 917, eq. 3).

function Mooncake.rrule!!(::CoDual{typeof(wrightomega)}, x::CoDual{Float64})
    w = wrightomega(primal(x))
    dw = w / (1.0 + w)
    wrightomega_pb!!(dy::Float64) = NoRData(), dw * dy
    return zero_fcodual(w), wrightomega_pb!!
end

function Mooncake.frule!!(::Dual{typeof(wrightomega)}, x::Dual{Float64})
    w = wrightomega(primal(x))
    return Dual(w, (w / (1.0 + w)) * tangent(x))
end

end
