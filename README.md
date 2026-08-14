# WrightOmega.jl

[![Build Status](https://github.com/jClugstor/WrightOmega.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jClugstor/WrightOmega.jl/actions/workflows/CI.yml?query=branch%3Amain)

The Wright omega function `ω(z)`: the solution of

```
ω + log(ω) = z
```

equal to `W_K(z)(exp(z))`, where `W` is the Lambert W function and `K(z)` is the
unwinding number. One exported function, `wrightomega`, for real and complex arguments.

## Usage

```julia
julia> using WrightOmega

julia> wrightomega(0.0)          # the omega constant
0.5671432904097838

julia> wrightomega(1.0)          # 1 + log(1) == 1
1.0

julia> wrightomega(1.0 + 1.0im)
0.9372082083733697 + 0.5054213160131512im

julia> wrightomega(1.0 + exp(1)) ≈ exp(1)       # 1.0 + exp(1) == log(exp(1)) + exp(1)
true
```



## Why not `lambertw(exp(x))`?

That composition overflows for `x > 709.78` and underflows for `x < -745.1`, even though
`ω(x)` itself is perfectly representable there (`ω(x) → x - log(x)` as `x → ∞`).
`wrightomega` never forms `exp(x)`, so it is finite and accurate over the entire
floating-point range:

```julia
julia> wrightomega(710.0)        # lambertw(exp(710.0)) would overflow
703.4440117119547

julia> wrightomega(1e12)
9.99999999972369e11
```

## Real arguments

Real evaluation follows Fukushima (2020): a piecewise minimax rational approximation on
each of twelve regions of the axis, with a rational in `exp(x)` below and one Newton
correction to `ω → x` above. The double-precision coefficients used here have a
measured maximum relative error of 8 machine epsilons.
The kernel flushes subnormal results, so `wrightomega(x) == 0.0` exactly for
`x < -708.4`.

## Complex arguments

Complex evaluation follows TOMS Algorithm 917 (Lawrence, Corless & Jeffrey 2012):
series initial approximations by region, polished by at most two fourth-order
iterations, accurate to a few ulp over the whole plane.

`ω` is single valued but discontinuous across the two rays `z = t ± iπ, t ≤ -1` (the
images of the Lambert W branch cut). Following the reference algorithm, an input whose
imaginary part is bitwise `±Float64(π)` with `real(z) ≤ -1` is taken to lie *on* the
ray, and the value there is the limit from below:

```julia
julia> wrightomega(-1 + im*π)                        # singular point
-1.0 + 0.0im

julia> wrightomega(complex(-2 + log(2), -Float64(π)))  # W₋₁ branch, below the cut
-2.0 - 0.0im
```

Inputs with zero imaginary part are dispatched to the real method, so
`wrightomega(x + 0.0im) == wrightomega(x)` exactly. The two code paths are different
algorithms; the test suite verifies they agree at and near the real axis.

## Derivatives

`dω/dz = ω/(1 + ω)`, and package extensions provide that rule to ForwardDiff, Enzyme,
Mooncake, ChainRulesCore (hence Zygote and friends), and register `wrightomega` with
Symbolics. Loading the corresponding package activates the extension; nothing else is
required.

## References

- T. Fukushima (2020). *Fast computation of the Wright ω function by piecewise minimax
  rational function approximation*. Real-argument method and coefficients.
- P. W. Lawrence, R. M. Corless, D. J. Jeffrey (2012). *Algorithm 917: Complex
  Double-Precision Evaluation of the Wright ω Function*. ACM Trans. Math. Software
  38(3). Complex-argument method.
- R. M. Corless, D. J. Jeffrey (2002). *The Wright ω Function*. Definition and
  properties.
