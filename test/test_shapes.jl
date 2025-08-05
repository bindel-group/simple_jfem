include("../src/shapes.jl")

using Test
using LinearAlgebra

function test_shape_lagrange(shapes)
    nen = nshapes(shapes)
    nodes = refnodes(shapes)
    L = zeros(nen, nen)
    for j = 1:nen
        L[:,j] .= shapes(nodes[:,j]).N[:]
    end
    @test L ≈ I
end

function test_dshape(shapes, x0)
    h = 1e-6
    for j = 1:dshapes(shapes)
        xp, xm = copy(x0), copy(x0)
        xp[j] = x0[j]+h
        xm[j] = x0[j]-h
        Np = copy(shapes(xp).N)
        Nm = copy(shapes(xm).N)
        @test shapes(x0).dN[:,j] ≈ (Np-Nm)/(2h)  atol=1e-8
    end
end

@testset "Lagrange property" begin
    test_shape_lagrange(Shapes1dP1())
    test_shape_lagrange(Shapes1dP2())
    test_shape_lagrange(Shapes1dP3())
    test_shape_lagrange(Shapes2dP1())
    test_shape_lagrange(Shapes2dP2())
    test_shape_lagrange(Shapes2dS2())
    test_shape_lagrange(Shapes2dT1())
end

@testset "1D shape derivatives" begin
    x1test = [0.876]
    test_dshape(Shapes1dP1(), x1test)
    test_dshape(Shapes1dP2(), x1test)
    test_dshape(Shapes1dP3(), x1test)
end

@testset "2D shape derivatives" begin
    x2test = [0.1312; 0.2488]
    test_dshape(Shapes2dP1(), x2test)
    test_dshape(Shapes2dP2(), x2test)
    test_dshape(Shapes2dS2(), x2test)
    test_dshape(Shapes2dT1(), x2test)
end
