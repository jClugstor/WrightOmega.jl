module WrightOmegaEnzymeExt

using Enzyme
using WrightOmega: wrightomega

# d/dx wrightomega(x) = Ω/(1 + Ω) where Ω = wrightomega(x)  (TOMS 917, eq. 3).
# @easy_rule binds Ω to the primal result and generates both the forward- and
# reverse-mode Enzyme rules; the constraint keeps the rule to real scalars.
Enzyme.EnzymeRules.@easy_rule(wrightomega(x::AbstractFloat),
    (Ω / (1 + Ω),)
)

end
