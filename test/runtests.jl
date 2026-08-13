using WrightOmega
using Test
using Symbolics: Symbolics, @variables
using ChainRulesCore: ChainRulesCore
using ForwardDiff: ForwardDiff
using Mooncake: Mooncake
using Enzyme: Enzyme
using Random: Xoshiro

const WO_PI = 3.141592653589793   # Float64(pi), bitwise the on-ray value

# Independent reference: refine the defining equation w + log(w) == z by Newton's
# method in 256-bit arithmetic, seeded from the double-precision result (so the branch
# is inherited; branch *selection* is tested separately via known values and the
# discontinuity tests below).
function refine(z::ComplexF64, w0::ComplexF64)
    setprecision(BigFloat, 256) do
        w = Complex{BigFloat}(w0)
        zb = Complex{BigFloat}(z)
        if real(w0) < 0 && abs(imag(w0)) < 1e-3 * abs(real(w0))
            # w sits (nearly) on log's branch cut, where plain Newton flips branches on
            # rounding noise. Iterate the regularized form -lam + log(lam) = z -+ i*pi
            # (lam = -w) instead, with the branch inherited from the seed's side.
            sigma = signbit(imag(w0)) ? -1 : 1
            zeta = zb - sigma * im * BigFloat(pi)
            lam = -w
            for _ in 1:8
                lam -= (-lam + log(lam) - zeta) * lam / (1 - lam)
            end
            return -lam
        end
        for _ in 1:8
            w -= (w + log(w) - zb) * w / (w + 1)
        end
        return w
    end
end

function relerr(z::ComplexF64)
    w = wrightomega(z)
    wref = refine(z, w)
    return Float64(abs(Complex{BigFloat}(w) - wref) / abs(wref))
end

@testset "WrightOmega.jl" begin

    @testset "real: values" begin
        @test wrightomega(0.0) ≈ 0.5671432904097838 rtol = 4eps()   # omega constant
        @test wrightomega(1.0) ≈ 1.0 rtol = 4eps()
        @test wrightomega(1.0 + ℯ) ≈ float(ℯ) rtol = 4eps()
        # defining equation over the normal range
        # (atol floor: w + log(w) cancels for small |x|, so the residual there reflects
        # rounding in the test expression, not error in w)
        for x in exp10.(range(-3, 300; length = 40))
            w = wrightomega(x)
            @test w + log(w) ≈ x rtol = 1e-13 atol = 1e-14
            w = wrightomega(-x)
            if w > 0                                 # not flushed to zero
                @test w + log(w) ≈ -x rtol = 1e-13 atol = 1e-14
            end
        end
        # monotonic on a coarse grid
        xs = -800.0:7.3:800.0
        @test issorted(wrightomega.(xs))
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
    end

    @testset "complex: special values (TOMS 917 Table I)" begin
        @test wrightomega(complex(-1.0, WO_PI)) ≈ -1.0 + 0.0im atol = 1e-15
        @test wrightomega(complex(-1.0, -WO_PI)) ≈ -1.0 + 0.0im atol = 1e-15
        # the two sides of the branch cut: W0 above, W_{-1} below
        @test wrightomega(complex(-2 + log(2), WO_PI)) ≈ -0.40637573995996 + 0im rtol = 1e-13
        @test wrightomega(complex(-2 + log(2), -WO_PI)) ≈ -2.0 + 0im rtol = 1e-13
        @test wrightomega(complex(0.0, 1 + WO_PI / 2)) ≈ im atol = 1e-15
        @test wrightomega(complex(0.0, WO_PI)) ≈
              -0.31813150520476413 + 1.3372357014306895im rtol = 1e-14   # W0(-1)
        @test wrightomega(1.0 + 1.0im) ≈ 0.9372082083733697 + 0.5054213160131512im rtol = 1e-14
    end

    @testset "complex: accuracy across all regions" begin
        pts = ComplexF64[
            -1.5+2.5im, 0.5+3im, 0.9+1.1im,             # region 1 (branch-point series)
            -1.5-2.5im, 0.5-3im,                        # region 2
            -3+1im, -10+3im, -300+0.5im, -2.1-3im,      # region 3 (series about -inf)
            0.5+0.5im, 1+1im, 2+2im, -1.99+0.3im,       # region 4 (Taylor at 1)
            1+1e-8im, 1.005+0.005im,                    # region 4, near-1 shortcut
            -5+3.5im, -100+10im, -2.5+3.2im,            # region 5 (negative-log series)
            -5-3.5im, -100-10im,                        # region 6
            10+10im, -5+8im, 100-50im, 3+7im,           # region 7 (asymptotic series)
            1e5+1e5im, 1e200+1e200im, 1e300-1e300im,    # region 7, huge |z|
            -5+(WO_PI-0.005)im, -5+(WO_PI+0.005)im,     # regularization band, top
            -5-(WO_PI-0.005)im, -5-(WO_PI+0.005)im,     # regularization band, bottom
            -0.995+WO_PI*im, -0.995-WO_PI*im,           # band, past the ray endpoints
            complex(-5.0, WO_PI), complex(-5.0, -WO_PI) # on the rays
        ]
        for z in pts
            @test relerr(z) < 20eps() # paper reports <= 2 ulp
        end
        # near the singular points -1 +- i*pi the condition number blows up like
        # 1/sqrt(distance); allow the correspondingly conditioned error
        for z in ComplexF64[-1+1e-6im+WO_PI*im, -1.000001-WO_PI*im]
            @test relerr(z) < 1e-12
        end
    end

    @testset "complex: functional identity w*exp(w) == exp(z)" begin
        for x in -17.0:1.7:17.0, y in -17.0:1.3:17.0
            abs(abs(y) - WO_PI) < 0.05 && continue   # skip the ill-conditioned rays
            z = complex(x, y)
            w = wrightomega(z)
            @test w * exp(w) ≈ exp(z) rtol = 1e-13
        end
    end

    @testset "complex: defining equation in the principal strip" begin
        for x in -30.0:3.1:30.0, y in (-WO_PI+0.05):0.5:WO_PI
            z = complex(x, y)
            w = wrightomega(z)
            @test w + log(w) ≈ z rtol = 1e-13
        end
    end

    @testset "complex: conjugate symmetry" begin
        for z in ComplexF64[1.5+2.3im, -3+1im, -5+3.5im, -1.5+2.5im, 10+10im,
                            -5+(WO_PI-0.003)im, 0.3+0.2im, -7.0+0.0im]
            @test wrightomega(conj(z)) == conj(wrightomega(z))
        end
    end

    @testset "complex: discontinuity across the rays" begin
        for t in (-2.0, -5.0, -20.0)
            # on the ray == limit from below the ray (continuity from below)
            w_on = wrightomega(complex(t, WO_PI))
            w_below = wrightomega(complex(t, prevfloat(WO_PI)))
            @test w_on ≈ w_below rtol = 1e-10
            @test real(w_on) < 0 && abs(imag(w_on)) < 1e-15    # W0-like small negative
            # crossing the ray jumps to the W_{-1}-like branch
            w_above = wrightomega(complex(t, nextfloat(WO_PI)))
            @test abs(w_above - w_on) > 1
            @test real(w_above) < -1 && imag(w_above) >= 0     # Im w > 0 above the ray
            # mirrored bottom ray: on-ray == limit from below (outside the strip)
            v_on = wrightomega(complex(t, -WO_PI))
            v_below = wrightomega(complex(t, prevfloat(-WO_PI)))
            @test v_on ≈ v_below rtol = 1e-10
            @test real(v_on) < -1 && abs(imag(v_on)) < 1e-15   # W_{-1}-like
            v_above = wrightomega(complex(t, nextfloat(-WO_PI)))
            @test abs(v_above - v_on) > 1
            @test imag(v_above) <= 0                           # Im w <= 0 in the strip
        end
    end

    @testset "complex: real-axis fast path routing" begin
        for x in (-100.0, -5.0, -1.5, 0.0, 0.5, 2.0, 1e6)
            @test wrightomega(complex(x, 0.0)) === complex(wrightomega(x), 0.0)
            @test wrightomega(complex(x, -0.0)) === complex(wrightomega(x), -0.0)
        end
    end

    @testset "complex: non-finite inputs" begin
        @test isnan(real(wrightomega(complex(NaN, 1.0))))
        @test isnan(imag(wrightomega(complex(1.0, NaN))))
        @test wrightomega(complex(-Inf, 0.0)) == 0.0 + 0.0im
        @test wrightomega(complex(-Inf, -0.5)) == 0.0 - 0.0im
        @test wrightomega(complex(-Inf, 2.0)) == 0.0 + 0.0im       # still inside the strip
        @test wrightomega(complex(-Inf, 7.0)) == complex(-Inf, 7.0 - WO_PI)
        @test wrightomega(complex(-Inf, -7.0)) == complex(-Inf, -7.0 + WO_PI)
        @test wrightomega(complex(Inf, 2.0)) == complex(Inf, 2.0)
        @test wrightomega(complex(2.0, Inf)) == complex(2.0, Inf)
    end

    @testset "complex: generic types" begin
        @test wrightomega(ComplexF32(1 + 1im)) isa ComplexF32
        @test wrightomega(ComplexF32(1 + 1im)) ≈ ComplexF32(wrightomega(1.0 + 1.0im))
        @test wrightomega(complex(1, 1)) === wrightomega(1.0 + 1.0im)   # Complex{Int}
        wb = wrightomega(Complex{BigFloat}(1 + 1im))                    # no recursion
        @test wb isa Complex{BigFloat}
        @test wb ≈ wrightomega(1.0 + 1.0im) rtol = 1e-14
    end

    @testset "Symbolics extension" begin
        @variables x
        expr = wrightomega(x)             # would MethodError if not registered
        @test expr isa Symbolics.Num
        x0 = 0.7
        w0 = wrightomega(x0)
        # compile to numeric functions (substitute-based constant folding differs
        # between Symbolics v6 and v7; build_function is stable across both)
        fw = Symbolics.build_function(expr, x; expression = Val(false))
        @test fw(x0) ≈ w0

        d = Symbolics.derivative(expr, x)
        fd = Symbolics.build_function(d, x; expression = Val(false))
        @test fd(x0) ≈ w0 / (1 + w0)
        # and against a central finite difference
        h = 1e-6
        cd = (wrightomega(x0 + h) - wrightomega(x0 - h)) / (2h)
        @test fd(x0) ≈ cd rtol = 1e-8
    end

    # The AD rules all encode the same identity, dω/dx = ω/(1+ω), computed through the
    # primal — so where a backend uses the rule, its result is bitwise equal to it.

    @testset "ChainRulesCore extension" begin
        x0 = 0.7
        w = wrightomega(x0)
        dw = w / (1 + w)
        Ω, pb = ChainRulesCore.rrule(wrightomega, x0)
        @test Ω == w
        f̄, x̄ = pb(1.0)
        @test f̄ === ChainRulesCore.NoTangent()
        @test x̄ == dw
        Ω2, Ω̇ = ChainRulesCore.frule((ChainRulesCore.NoTangent(), 1.0), wrightomega, x0)
        @test Ω2 == w
        @test Ω̇ == dw
    end

    @testset "ForwardDiff extension" begin
        x0 = 0.7
        w = wrightomega(x0)
        @test ForwardDiff.derivative(wrightomega, x0) == w / (1 + w)
        # second derivative through nested Duals; analytically w/(1+w)^3
        d2 = ForwardDiff.derivative(x -> ForwardDiff.derivative(wrightomega, x), x0)
        @test d2 ≈ w / (1 + w)^3 rtol = 1e-14
        @test ForwardDiff.derivative(wrightomega, 0.7f0) isa Float32
        @test ForwardDiff.derivative(wrightomega, 0.7f0) ≈ Float32(w / (1 + w))
    end

    @testset "Mooncake extension" begin
        rng = Xoshiro(123)
        # correctness vs finite differences, both forward and reverse mode
        for x in (0.7, -5.0, 40.0)
            Mooncake.TestUtils.test_rule(rng, wrightomega, x; print_results = false)
        end
        cache = Mooncake.prepare_gradient_cache(wrightomega, 0.7)
        v, g = Mooncake.value_and_gradient!!(cache, wrightomega, 0.7)
        w = wrightomega(0.7)
        @test v == w
        @test g[2] == w / (1 + w)
    end

    @testset "Enzyme extension" begin
        x0 = 0.7
        w = wrightomega(x0)
        dw = w / (1 + w)
        @test Enzyme.autodiff(Enzyme.Forward, wrightomega, Enzyme.Duplicated(x0, 1.0))[1] == dw
        @test Enzyme.autodiff(Enzyme.Reverse, wrightomega, Enzyme.Active(x0))[1][1] == dw
    end

end
