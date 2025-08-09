using DispatchDoctor: @stable

@stable begin
include("../src/quadrules.jl")
end

using Test

@testset "Test 1D quadrature" begin
    ptest(xi) = ((((xi+3)*xi-1)*xi+1)*xi+1)*xi+1
    Iref = 6.0/5.0 + 2.0/3.0 + 2.0
    for npts = 3:10
        I = sum(ptest(xi) * wt for (xi, wt) in GaussRule1d(npts))
        @test I ≈ Iref  atol=1e-12
    end
end

@testset "Test 2D Gauss" begin
    for d = 1:5
        Iref = (4.0, 0.0, 0.0, 0.0)
        for pxy = 1:4
            px = (pxy-1)%2
            py = div(pxy-1,2)
            f(xi) = xi[1]^px * xi[2]^py
            I = sum(f(xi) * wt for (xi, wt) in GaussRule2d(d))
            @test I ≈ Iref[pxy]  atol=1e-12
        end
    end
end

@testset "Hughes rule for triangles" begin
    Iref = (0.5, 1.0/6, 1.0/6, 1.0/24)
    for pxy = 1:4
        px = (pxy-1)%2
        py = div(pxy-1,2)
        f(xi) = xi[1]^px * xi[2]^py
        I = sum(f(xi) * wt for (xi, wt) in HughesRule2d())
        @test I ≈ Iref[pxy]  atol=1e-12
    end
end
