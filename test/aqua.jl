@testset "Aqua" begin
    # ambiguities/piracy also cover the extension methods, since the AD packages
    # are loaded in runtests.jl
    Aqua.test_all(WrightOmega)
end
