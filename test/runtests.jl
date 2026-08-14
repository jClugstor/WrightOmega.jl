using WrightOmega
using Test
using Symbolics: Symbolics, @variables
using ChainRulesCore: ChainRulesCore
using ForwardDiff: ForwardDiff
using Mooncake: Mooncake
using Enzyme: Enzyme
using Random: Xoshiro
using Aqua: Aqua
using JET: JET

const WO_PI = 3.141592653589793   # Float64(pi), bitwise the on-ray value

@testset "WrightOmega.jl" begin
    include("aqua.jl")
    include("jet.jl")
    include("real.jl")
    include("complex.jl")
    include("extensions.jl")
end
