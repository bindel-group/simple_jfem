include("../src/shapes.jl")

using Test
using LinearAlgebra

function test_shape_lagrange(shape!, nodes)
    d, nen = size(nodes)
    N  = zeros(nen)
    dN = zeros(nen, d)
    L  = zeros(nen, nen)
    for j = 1:nen
        shape!(N, dN, nodes[:,j])
        L[:,j] .= N[:]
    end
    @test L ≈ I
end

function test_dshape(shape!, x0, nen)
    d = length(x0)
    N, Np, Nm = zeros(nen), zeros(nen), zeros(nen)
    dN = zeros(nen, d)
    h  = 1e-6
    for j = 1:d
        xp, xm = copy(x0), copy(x0)
        xp[j] = x0[j]+h
        xm[j] = x0[j]-h
        shape!(Np, dN, xp)
        shape!(Nm, dN, xm)
        shape!(N, dN, x0)
        @test dN[:,j] ≈ (Np-Nm)/(2h)  atol=1e-8
    end
end

@testset "Lagrange property" begin
test_shape_lagrange(shapes1dP1!, [-1.0 1.0])
test_shape_lagrange(shapes1dP2!, [-1.0 0.0 1.0])
test_shape_lagrange(shapes1dP3!, [-1.0 -1.0/3 1.0/3 1.0])
test_shape_lagrange(shapes2dP1!,
        [-1.0   1.0   1.0  -1.0 ;
         -1.0  -1.0   1.0   1.0])
test_shape_lagrange(shapes2dP2!,
        [-1.0   0.0   1.0   1.0   1.0   0.0  -1.0  -1.0   0.0 ;
         -1.0  -1.0  -1.0   0.0   1.0   1.0   1.0   0.0   0.0 ])
test_shape_lagrange(shapes2dS2!,
        [-1.0   0.0   1.0   1.0   1.0   0.0  -1.0  -1.0 ;
         -1.0  -1.0  -1.0   0.0   1.0   1.0   1.0   0.0 ])
end

@testset "1D shape derivatives" begin
x1test = [0.876]
test_dshape(shapes1dP1!, x1test, 2)
test_dshape(shapes1dP2!, x1test, 3)
test_dshape(shapes1dP3!, x1test, 4)
end

@testset "2D shape derivatives" begin
x2test = [0.1312; 0.2488]
test_dshape(shapes2dP1!, x2test, 4)
test_dshape(shapes2dP2!, x2test, 9)
test_dshape(shapes2dS2!, x2test, 8)
test_dshape(shapes2dT1!, x2test, 3)
end
