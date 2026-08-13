# Complex evaluation of the Wright omega function.
#
# Method: Lawrence, Corless & Jeffrey (2012), "Algorithm 917: Complex Double-Precision
# Evaluation of the Wright omega Function", ACM TOMS 38(3). A series initial
# approximation chosen by region of the complex plane (their Algorithm 7.1), followed by
# at most two iterations of the fourth-order FSC scheme, eq. (15) with N = 3.
#
# omega is single valued but discontinuous across the rays z = t +- i*pi, t <= -1 (the
# images under z -> exp(z) of the Lambert W branch cut). Following the paper, an input
# whose imaginary part is bitwise +-Float64(pi) with real part <= -1 is declared to lie
# ON a ray, where the function is continuous from below. Near the rays the defining
# equation omega + log(omega) = z is regularized to -lambda + log(lambda) = z -+ i*pi
# with lambda = -omega (paper section 6), so the iteration never evaluates log across
# its branch cut. The paper prescribes directed rounding for the z -+ i*pi subtraction;
# none is needed here: inside the |Im(z) -+ pi| <= 0.01 band the subtraction is exact by
# Sterbenz's lemma, and the only ambiguity -- Im(z) bitwise equal to +-pi, which must
# behave as approached from below -- is resolved by flipping the difference +0.0 to -0.0.

# The kernel must compare against the *rounded* pi. Comparisons against the Irrational
# see the true value, so e.g. `y > -pi` would be true for y == -Float64(pi), silently
# putting on-ray inputs in the wrong region.
const WO_PI = 3.141592653589793

# Negative-log series, paper eq. (27), used just outside the rays with t = z -+ i*pi.
# log(-t) is continuous there (Re(-t) >= 2), which is the point of this form.
@inline function _wo_neglog(t::ComplexF64)
    t1  = t + 1.0
    L   = log(-t)
    t12 = t1 * t1
    t15 = t12 * t12 * t1
    return t - t * L / t1 + t * L^2 / (2.0 * t12 * t1) +
           t * (2.0 * t - 1.0) * L^3 / (6.0 * t15) +
           t * (6.0 * t^2 - 8.0 * t + 1.0) * L^4 / (24.0 * t15 * t12) +
           t * (24.0 * t^3 - 58.0 * t^2 + 22.0 * t - 1.0) * L^5 / (120.0 * t15 * t12 * t12)
end

# One FSC step (eq. 15 with N = 3, rearranged as in section 7.1 of the paper). The
# correction ratio (q - r/2)/(q - r) is written with q divided out so that when q
# overflows (|w| beyond ~1e154) the ratio tends to 1 and the step degrades gracefully
# to a Newton step instead of producing Inf/Inf.
@inline function _wo_fsc_step(w::ComplexF64, w1::ComplexF64, r::ComplexF64)
    q = w1 * (w1 + (2.0 / 3.0) * r)
    ratio = (1.0 - r / (2.0 * q)) / (1.0 - r / q)
    # |w1| beyond ~1e154 overflows q (to NaN via Inf - Inf when Re and Im are of the
    # same magnitude); the true ratio is then 1 to full precision.
    isfinite(real(ratio)) && isfinite(imag(ratio)) || (ratio = complex(1.0, 0.0))
    return w * (1.0 + (r / w1) * ratio)
end

function _wrightomega(z::ComplexF64)
    x, y = reim(z)

    # Signed distance to each ray; exact by Sterbenz inside the regularization band.
    ympi = y - WO_PI
    yppi = y + WO_PI
    # On-ray inputs (bitwise +-pi, x <= -1) take the limit from below the ray.
    if x <= -1.0
        y == WO_PI && (ympi = -0.0)
        y == -WO_PI && (yppi = -0.0)
    end

    # ------------------------------------------------------ initial approximation
    local w::ComplexF64
    if -2.0 < x <= 1.0 && 1.0 < y < 2.0 * WO_PI
        # Region 1: branch-point series about -1 + i*pi, eq. (21). The imaginary part
        # of 2(z + 1 - i*pi) is built from ympi so the on-ray sign flip above selects
        # the correct square root on the sqrt branch cut.
        p = sqrt(complex(2.0 * (x + 1.0), 2.0 * ympi))
        w = -1.0 + p * (im + p * (1 / 3 + p * (-im / 36 + p * (1 / 270 + p * (im / 4320)))))
    elseif -2.0 < x <= 1.0 && -2.0 * WO_PI < y < -1.0
        # Region 2: conjugate branch-point series about -1 - i*pi, eq. (23).
        p = sqrt(complex(2.0 * (x + 1.0), 2.0 * yppi))
        w = -1.0 + p * (-im + p * (1 / 3 + p * (im / 36 + p * (1 / 270 + p * (-im / 4320)))))
    elseif x <= -2.0 && -WO_PI < y <= WO_PI
        # Region 3: series about -inf inside the strip, eq. (24), in powers of exp(z).
        t = exp(z)
        w = t * (1.0 + t * (-1.0 + t * (1.5 + t * (-8 / 3 + t * (125 / 24)))))
    elseif (-2.0 < x <= 1.0 && -1.0 <= y <= 1.0) ||
           (x > -1.0 && (x - 1.0)^2 + y^2 <= WO_PI^2)
        # Region 4: Taylor series about the regular point z = 1, eq. (29).
        u = z - 1.0
        w = 1.0 + u * (1 / 2 + u * (1 / 16 + u * (-1 / 192 + u * (-1 / 3072 + u * (13 / 61440)))))
        # Truly near z = 1 the series alone is already at full precision and the
        # residual z - w - log(w) is ill-conditioned (w ~ 1); return it directly
        # (paper section 4.5).
        abs2(u) < 1.0e-4 && return w
    elseif x <= -2.0 && WO_PI < y && ympi <= -0.75 * (x + 1.0)
        # Region 5: negative-log series just above the top ray, t = z - i*pi.
        w = _wo_neglog(complex(x, ympi))
    elseif x <= -2.0 && y <= -WO_PI && yppi >= 0.75 * (x + 1.0)
        # Region 6: mirror of region 5, just below the bottom ray, t = z + i*pi.
        w = _wo_neglog(complex(x, yppi))
    else
        # Region 7: asymptotic series about infinity, eq. (25). Powers of z are built
        # by repeated (overflow-robust) division; z*z itself can overflow to NaN.
        L = log(z)
        Lz = L / z
        Lz2 = Lz / z
        Lz3 = Lz2 / z
        w = z - L + Lz + Lz2 * (0.5 * L - 1.0) + Lz3 * (L * L / 3 - 1.5 * L + 1.0)
    end

    # ------------------------------------------------------ regularization (section 6)
    # Near the rays, solve -lambda + log(lambda) = z -+ i*pi for lambda = -omega instead.
    # Since log(s*w) and log(w) differ locally by a constant, the FSC update is
    # unchanged; only the residual differs: r = zeta - w - log(s*w) with s = -1 and
    # shifted zeta. lambda = s*w stays near the positive real axis, where log is
    # continuous, so rounding can never kick the iteration across the discontinuity.
    s = 1.0
    zeta = z
    if x <= -0.99 && (abs(ympi) <= 0.01 || abs(yppi) <= 0.01)
        s = -1.0
        zeta = complex(x, abs(ympi) <= 0.01 ? ympi : yppi)
    end

    # ------------------------------------------------------ FSC iteration (at most two)
    # exp(z) underflowed between the rays (x < ~-745); omega = 0 is the correct limit.
    iszero(w) && return w
    r = zeta - w - log(s * w)
    w1 = w + 1.0
    # r == 0: the series is already exact (in particular at the singular points
    # -1 +- i*pi, where w == -1 and w1 == 0 would otherwise poison the update).
    if !iszero(r) && !iszero(w1)
        w = _wo_fsc_step(w, w1, r)
        # Second step only if the predicted next residual, eq. (19), is above eps.
        if abs(r)^4 * abs(2.0 * w * w - 8.0 * w - 1.0) >= 72.0 * eps(Float64) * abs(w1)^6
            r = zeta - w - log(s * w)
            w1 = w + 1.0
            if !iszero(r) && !iszero(w1)
                w = _wo_fsc_step(w, w1, r)
            end
        end
    end

    # On the rays omega is exactly real (and negative), so any residual imaginary part
    # is rounding noise; more importantly its zero must be signed to match the
    # limit-from-below convention, so that w + log(w) == z holds with the principal
    # log: +0.0 on the top ray (log branch +i*pi), -0.0 on the bottom (-i*pi).
    if x <= -1.0
        y == WO_PI && return complex(real(w), 0.0)
        y == -WO_PI && return complex(real(w), -0.0)
    end
    return w
end

"""
    wrightomega(z::Complex) -> Complex

The Wright omega function on the complex plane: `omega(z) = W_K(z)(exp(z))` with `K` the
unwinding number -- equivalently, the solution of `w + log(w) == z` continued from the
real value `omega(0) = 0.567...`.

Single valued everywhere, but discontinuous across the two rays `z = t ± im*π, t ≤ -1`
(the images of the Lambert `W` branch cut). Following Lawrence, Corless & Jeffrey
(2012), an input whose imaginary part equals `±Float64(π)` bitwise with `real(z) ≤ -1`
lies *on* a ray, and the value there is the limit from below the ray. Inputs with zero
imaginary part are real and take the fast real path.

Accurate to a few ulp over the whole plane: a series initial approximation by region is
polished by at most two fourth-order iterations of the defining equation.

```
julia> wrightomega(1.0 + 0.0im)
1.0 + 0.0im

julia> wrightomega(-1 + im*π)      # singular point, omega == -1
-1.0 + 0.0im
```

Method: TOMS Algorithm 917 (Lawrence, Corless & Jeffrey 2012). See `src/complex.jl`.
"""
function wrightomega(z::ComplexF64)
    x, y = reim(z)
    if isnan(x) || isnan(y)
        return ComplexF64(NaN, NaN)
    elseif isinf(x) && x < 0.0 && -WO_PI < y <= WO_PI
        # omega -> 0 as Re(z) -> -inf inside the strip between the rays.
        return complex(0.0, copysign(0.0, y))
    elseif isinf(x) && x < 0.0 && isfinite(y)
        # Outside the strip omega ~ z -+ i*pi (the region-5/6 asymptote).
        return complex(x, y - copysign(WO_PI, y))
    elseif isinf(x) || isinf(y)
        # omega ~ z toward every other infinity.
        return z
    elseif iszero(y)
        # omega is real on the whole real axis: take the fast real path. The zero
        # imaginary part keeps its sign so that wrightomega(conj(z)) == conj(wrightomega(z)).
        return complex(wrightomega(x), y)
    else
        return _wrightomega(z)
    end
end

wrightomega(z::T) where {T<:Complex} = oftype(float(z), wrightomega(ComplexF64(z)))
