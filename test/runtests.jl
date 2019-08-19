using Test

using HarmonieVerticalLevels

@testset "A B" begin
    @test A(0.0) == 0.0
    @test A(1.0) == 0.0
    @test B(0.0) == 0.0
    @test B(1.0) == 1.0  
end

