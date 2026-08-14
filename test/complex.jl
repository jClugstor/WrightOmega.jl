# Tests for the complex path (TOMS Algorithm 917) and its seam with the real path.

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
    @test all(relerr(z) < 20eps() for z in pts)   # paper reports <= 2 ulp
    # near the singular points -1 +- i*pi the condition number blows up like
    # 1/sqrt(distance); allow the correspondingly conditioned error
    @test all(relerr(z) < 1e-12 for z in ComplexF64[-1+1e-6im+WO_PI*im, -1.000001-WO_PI*im])
end

@testset "complex: functional identity w*exp(w) == exp(z)" begin
    function identity_ok(z)
        abs(abs(imag(z)) - WO_PI) < 0.05 && return true   # skip the ill-conditioned rays
        w = wrightomega(z)
        return isapprox(w * exp(w), exp(z); rtol = 1e-13)
    end
    @test all(identity_ok(complex(x, y)) for x in -17.0:1.7:17.0, y in -17.0:1.3:17.0)
end

@testset "complex: defining equation in the principal strip" begin
    function strip_ok(z)
        w = wrightomega(z)
        return isapprox(w + log(w), z; rtol = 1e-13)
    end
    @test all(strip_ok(complex(x, y))
              for x in -30.0:3.1:30.0, y in (-WO_PI+0.05):0.5:WO_PI)
end

@testset "complex: conjugate symmetry" begin
    @test all(wrightomega(conj(z)) == conj(wrightomega(z))
              for z in ComplexF64[1.5+2.3im, -3+1im, -5+3.5im, -1.5+2.5im, 10+10im,
                                  -5+(WO_PI-0.003)im, 0.3+0.2im, -7.0+0.0im])
end

@testset "complex: discontinuity across the rays" begin
    function rays_ok(t)
        # top ray: on-ray == limit from below (inside the strip, W0-like small
        # negative real value); crossing it jumps to the W_{-1}-like branch
        w_on = wrightomega(complex(t, WO_PI))
        w_below = wrightomega(complex(t, prevfloat(WO_PI)))
        w_above = wrightomega(complex(t, nextfloat(WO_PI)))
        top = isapprox(w_on, w_below; rtol = 1e-10) &&
              real(w_on) < 0 && abs(imag(w_on)) < 1e-15 &&
              abs(w_above - w_on) > 1 &&
              real(w_above) < -1 && imag(w_above) >= 0     # Im w > 0 above the ray
        # bottom ray, mirrored: on-ray == limit from below (outside the strip)
        v_on = wrightomega(complex(t, -WO_PI))
        v_below = wrightomega(complex(t, prevfloat(-WO_PI)))
        v_above = wrightomega(complex(t, nextfloat(-WO_PI)))
        bot = isapprox(v_on, v_below; rtol = 1e-10) &&
              real(v_on) < -1 && abs(imag(v_on)) < 1e-15 &&
              abs(v_above - v_on) > 1 &&
              imag(v_above) <= 0                           # Im w <= 0 in the strip
        return top && bot
    end
    @test all(rays_ok(t) for t in (-2.0, -5.0, -20.0))
end

@testset "complex: real-axis fast path routing" begin
    @test all(wrightomega(complex(x, y)) === complex(wrightomega(x), y)
              for x in (-100.0, -5.0, -1.5, 0.0, 0.5, 2.0, 1e6), y in (0.0, -0.0))
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

@testset "real/complex seam: near-axis consistency" begin
    # The real path (Fukushima) and the complex path (TOMS 917) are different
    # algorithms; a few ulp off the real axis the complex kernel must agree with
    # the real kernel to their combined accuracy. To first order
    # ω(x + iε) = ω(x) + iε ω/(1+ω) with 0 < ω/(1+ω) < 1, so the real part matches
    # ω(x) and the imaginary part is bounded by |ε| and carries ε's sign.
    function seam_consistent(x, s)
        ε = s * eps(x)
        wr = wrightomega(x)
        wc = wrightomega(complex(x, ε))
        return isapprox(real(wc), wr; rtol = 1e-14) &&
               abs(imag(wc)) <= 2 * abs(ε) &&
               (sign(imag(wc)) == sign(ε) || iszero(imag(wc)))
    end
    @test all(seam_consistent(x, s)
              for x in (-500.0, -5.0, -1.7969, -1.0, 0.0, 1.0, 2.7, 10.0, 40.0,
                        250.0, 1.5e3, 1e4, 2e5, 2e6, 1e8, 1e10, 4e11, 1e12),
                  s in (4.0, -4.0))
end
