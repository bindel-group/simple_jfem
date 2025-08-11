using DispatchDoctor: @stable

@stable begin
include("../src/shapes.jl")
include("../src/quadrules.jl")
include("../src/eltmap.jl")
end

using Test

@testset "Test mapping" begin
    s = Shapes2dP1()
    q = GaussRule2d(2)

    xyref = [0.5; 1.0]
    X = [ 0.0 1.0 1.0 0.0;
          0.0 0.0 1.0 1.0 ]
    J = zeros(2,2)
    ipiv = zeros(Int,2)

    # Checked mapped integration with explicit map
    function chi!(xy, J)
        x, y = xy
        xy[1] = 3.0 + 2*x
        xy[2] = 1.0 + x + y
        J[:,:] .= [2.0 0.0; 1.0 1.0]
    end
    I = sum( wt for (x,wt) in MappedRule(q, chi!) )
    @test I ≈ 8.0

    # Trivial mapping
    xy = isoparametric!(s, X, copy(xyref), J)
    r2s_factor!(J, ipiv)
    @test xy == [0.75; 1.0]
    @test r2s_det(J,ipiv) == 0.25

    # Check integration
    I = sum( wt for (x,wt) in IsoMappedRule(q, s, X) )
    @test I ≈ 1.0

    # More interesting mapping
    emap(x, y) = (3.0+2*x, 1.0+x+y)
    for i = 1:4
        X[:,i] .= emap(X[:,i]...)
    end
    xy = isoparametric!(s, X, copy(xyref), J)
    r2s_factor!(J, ipiv)
    @test xy == [4.5; 2.75]
    @test r2s_det(J, ipiv) == 0.5

    # Check integration
    I = sum( wt for (x,wt) in IsoMappedRule(q, s, X) )
    @test I ≈ 2.0
end
