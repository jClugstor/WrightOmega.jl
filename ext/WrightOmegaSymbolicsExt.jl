module WrightOmegaSymbolicsExt

# import (not using): @register_symbolic adds methods to wrightomega
import WrightOmega: wrightomega
using Symbolics: Symbolics, @register_symbolic

# Keep wrightomega(::Num) as a symbolic call node instead of tracing into the Float64
# kernel (whose generic Real method would error on symbolic input anyway).
@register_symbolic wrightomega(x)

# d/dx wrightomega(x) = w/(1 + w) with w = wrightomega(x)  (TOMS 917, eq. 3).
# Expressed through wrightomega itself so the derivative shares the primal evaluation.
@static if pkgversion(Symbolics) >= v"7"
    Symbolics.@register_derivative wrightomega(x) 1 wrightomega(x) / (1 + wrightomega(x))
else
    # Symbolics v6 predates @register_derivative; extend the dispatch hook directly.
    # wrap/unwrap: the rule receives unwrapped symbolics and must return one.
    function Symbolics.derivative(::typeof(wrightomega), args::NTuple{1,Any}, ::Val{1})
        w = wrightomega(Symbolics.wrap(args[1]))
        return Symbolics.unwrap(w / (1 + w))
    end
end

end
