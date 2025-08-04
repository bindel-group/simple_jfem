include("../src/quadrules.jl")

using Test

@testset "Test 1D quadrature" begin
    ptest(xi) = ((((xi+3)*xi-1)*xi+1)*xi+1)*xi+1
    Iref = 6.0/5.0 + 2.0/3.0 + 2.0

    for npts = 3:10
        I = 0.0
        for i = 1:npts
            I += ptest(gauss_point(i, npts)) * gauss_weight(i, npts)
        end
        @test I ≈ Iref  atol=1e-12
    end
end

@testset "Test 2D Gauss" begin
    for d = 1:5
        npts = d*d
        Iref = (4.0, 0.0, 0.0, 0.0)
        for pxy = 1:4
            px = (pxy-1)%2
            py = div(pxy-1,2)

            I = 0.0
            for i = 1:npts
                xi = gauss2d_point!(zeros(2), i, npts)
                wt = gauss2d_weight(i, npts)
                I += xi[1]^px * xi[2]^py * wt
            end

            @test I ≈ Iref[pxy]  atol=1e-12
        end
    end
end

@testset "Hughes rule for triangles" begin
    Iref = (0.5, 1.0/6, 1.0/6, 1.0/24)
    for pxy = 1:4
        px = (pxy-1)%2
        py = div(pxy-1,2)

        I = 0.0
        for i = 1:3
            xi = hughes_point!(zeros(2), i, 3)
            wt = hughes_weight(i, 3)
            I += xi[1]^px * xi[2]^py * wt
        end
        @test I ≈ Iref[pxy]  atol=1e-12
    end
end
