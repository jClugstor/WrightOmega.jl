# Tests for the package extensions (Symbolics + the AD backends).
#
# The AD rules all encode the same identity, dω/dx = ω/(1+ω), computed through the
# primal — so where a backend uses the rule, its result is bitwise equal to it.

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
