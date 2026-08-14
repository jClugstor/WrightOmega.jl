@testset "JET" begin
    # JET tracks compiler internals closely; run its checks only on released Julia
    # versions it targets, and let the other CI jobs cover functionality.
    if VERSION >= v"1.12" && isempty(VERSION.prerelease)
        # whole-package static analysis: no possible MethodErrors/undefined bindings
        JET.test_package(WrightOmega)
        # no runtime dispatch in the hot paths
        JET.@test_opt wrightomega(1.0)
        JET.@test_opt wrightomega(1.0 + 1.0im)
        JET.@test_opt wrightomega(complex(-5.0, 3.14))   # regularization band
    end
end
