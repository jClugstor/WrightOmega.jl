# Tests for the real path (Fukushima 2020 piecewise minimax rationals).

@testset "real: values" begin
    @test wrightomega(0.0) ≈ 0.5671432904097838 rtol = 4eps()   # omega constant
    @test wrightomega(1.0) ≈ 1.0 rtol = 4eps()
    @test wrightomega(1.0 + ℯ) ≈ float(ℯ) rtol = 4eps()
    # defining equation over the normal range
    # (atol floor: w + log(w) cancels for small |x|, so the residual there reflects
    # rounding in the test expression, not error in w; w == 0 is the documented
    # flush to zero for very negative x)
    function defeq(x)
        w = wrightomega(x)
        return w == 0 || isapprox(w + log(w), x; rtol = 1e-13, atol = 1e-14)
    end
    @test all(defeq(s * x) for x in exp10.(range(-3, 300; length = 1000)), s in (1, -1))
end

@testset "real: extremes and non-finite" begin
    @test wrightomega(1e300) ≈ 1e300 rtol = 1e-12       # ~ x - log(x)
    @test wrightomega(-1e300) == 0.0                    # documented flush to zero
    @test wrightomega(Inf) == Inf
    @test wrightomega(-Inf) == 0.0
    @test isnan(wrightomega(NaN))
end

@testset "real: generic types" begin
    @test wrightomega(1.0f0) === Float32(wrightomega(1.0))
    @test wrightomega(Float16(1.0)) isa Float16
    @test wrightomega(1) === wrightomega(1.0)
    @test wrightomega(1 // 2) === wrightomega(0.5)
    @test wrightomega(true) === wrightomega(1.0)
    wb = wrightomega(big"1.0")                          # no infinite recursion
    @test wb isa BigFloat && wb == 1.0
    @test wrightomega(missing) === missing              # like Base/SpecialFunctions
end
