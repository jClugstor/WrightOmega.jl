module WrightOmegaChainRulesCoreExt

using ChainRulesCore: @scalar_rule
using WrightOmega: wrightomega

# d/dx wrightomega(x) = Ω/(1 + Ω) where Ω = wrightomega(x)  (TOMS 917, eq. 3).
# @scalar_rule binds Ω to the primal output and derives both frule and rrule.
@scalar_rule(wrightomega(x), (Ω / (1 + Ω),))

end
