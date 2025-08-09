using DispatchDoctor: @stable

@stable begin
include("../src/shapes.jl")
include("../src/mesh.jl")
end

using Test

# Hand checked reference mesh for P1 quads

XrefP1 = [ 0.0    0.0 ;
           1.0/3  0.0 ;
           2.0/3  0.0 ;
           1.0    0.0 ;
           0.0    0.5 ;
           1.0/3  0.5 ;
           2.0/3  0.5 ;
           1.0    0.5 ;    
           0.0    1.0 ;
           1.0/3  1.0 ;
           2.0/3  1.0 ;
           1.0    1.0 ]'

erefP1 = [ 0 1 5  4  ;
           1 2 6  5  ;
           2 3 7  6  ;
           4 5 9  8  ;
           5 6 10 9  ;
           6 7 11 10 ]' .+ 1

# Hand checked reference mesh for P2 quads

XrefP2 = [ 0.0/6  0.0  ;
           1.0/6  0.0  ;
           2.0/6  0.0  ;
           3.0/6  0.0  ;
           4.0/6  0.0  ;
           5.0/6  0.0  ;
           1      0.0  ;
           0.0/6  0.25 ;
           1.0/6  0.25 ;
           2.0/6  0.25 ;
           3.0/6  0.25 ;
           4.0/6  0.25 ;
           5.0/6  0.25 ;
           1      0.25 ;
           0.0/6  0.5  ;
           1.0/6  0.5  ;
           2.0/6  0.5  ;
           3.0/6  0.5  ;
           4.0/6  0.5  ;
           5.0/6  0.5  ;
           1      0.5  ;
           0.0/6  0.75 ;
           1.0/6  0.75 ;
           2.0/6  0.75 ;
           3.0/6  0.75 ;
           4.0/6  0.75 ;
           5.0/6  0.75 ;
           1      0.75 ;
           0.0/6  1.0  ;
           1.0/6  1.0  ;
           2.0/6  1.0  ;
           3.0/6  1.0  ;
           4.0/6  1.0  ;
           5.0/6  1.0  ;
           1      1.0  ]'

erefP2 = [ 0    1    2    9   16   15   14    7    8  ;
           2    3    4   11   18   17   16    9   10  ;
           4    5    6   13   20   19   18   11   12  ;
          14   15   16   23   30   29   28   21   22  ;
          16   17   18   25   32   31   30   23   24  ;
          18   19   20   27   34   33   32   25   26  ]' .+ 1

# Hand checked reference mesh for S2 quads

XrefS2 = [ 0.0/6  0.0  ;
           1.0/6  0.0  ;
           2.0/6  0.0  ;
           3.0/6  0.0  ;
           4.0/6  0.0  ;
           5.0/6  0.0  ;
           1.0    0.0  ;
           0.0/3  0.25 ;
           1.0/3  0.25 ;
           2.0/3  0.25 ;
           1.0    0.25 ;
           0.0/6  0.5  ;
           1.0/6  0.5  ;
           2.0/6  0.5  ;
           3.0/6  0.5  ;
           4.0/6  0.5  ;
           5.0/6  0.5  ;
           1.0    0.5  ;
           0.0/3  0.75 ;
           1.0/3  0.75 ;
           2.0/3  0.75 ;
           1.0    0.75 ;
           0.0/6  1.0  ;
           1.0/6  1.0  ;
           2.0/6  1.0  ;
           3.0/6  1.0  ;
           4.0/6  1.0  ;
           5.0/6  1.0  ;
           1.0    1.0  ]'

erefS2 = [ 0   1   2   8  13  12  11   7 ;
           2   3   4   9  15  14  13   8 ;
           4   5   6  10  17  16  15   9 ;
          11  12  13  19  24  23  22  18 ;
          13  14  15  20  26  25  24  19 ;
          15  16  17  21  28  27  26  20 ]' .+ 1

# Hand checked reference mesh for triangles (nodes = P1 quad nodes)

erefT1 = [ 0   1   4 ;
           4   1   5 ;
           1   2   5 ;
           5   2   6 ;
           2   3   6 ;
           6   3   7 ;
           4   5   8 ;
           8   5   9 ;
           5   6   9 ;
           9   6  10 ;
           6   7  10 ;
          10   7  11 ]' .+ 1

function check_mesh(mesh, Xref, eref)
    @test mesh.X ≈ Xref atol=1e-8
    @test mesh.elt == eref
end

@testset "Check block meshers" begin
    check_mesh(mesh_block2d_P1(3,2), XrefP1, erefP1)
    check_mesh(mesh_block2d_P2(3,2), XrefP2, erefP2)
    check_mesh(mesh_block2d_S2(3,2), XrefS2, erefS2)
    #check_mech(mesh_block2d_T1(3,2), XrefP1, erefT1)
end

@testset "Test mapping" begin
    mesh = mesh_block2d_P1(1,1)
    xyref = [0.5; 1.0]

    # Trivial mapping
    xy, detJ = mesh_to_spatial(mesh, 1, xyref)
    @test xy == [0.75; 1.0]
    @test detJ == 0.25
    #@test F.factors == [0.5 0.0 ; 0.0 0.5]

    # More interesting mapping
    emap(x, y) = (3.0+2*x, 1.0+x+y)
    for i = 1:4
        mesh.X[:,i] .= emap(mesh.X[:,i]...)
    end
    xy, detJ = mesh_to_spatial(mesh, 1, xyref)
    @test xy == [4.5; 2.75]
    @test detJ == 0.5
    #@test F.factors == [1.0 0.0; 0.5 0.5]
end
