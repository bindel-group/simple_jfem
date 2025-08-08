include("../src/quadrules.jl")
include("../src/shapes.jl")
include("../src/mesh.jl")
include("../src/assemble.jl")
include("../src/fem.jl")
include("../src/element.jl")

using Test

function setup_test_mesh(numelt, shapes, qrule, u0, u1)
    mesh = mesh_create1d(numelt, shapes)
    fe = FEMProblem(mesh, Poisson1dElt(), qrule, 1)
    fe.id[1,1]   = -1
    fe.id[1,end] = -1
    fe.U[1,1]    = u0
    fe.U[1,end]  = u1
    assign_ids!(fe)
    fe
end

function test_fem1(shapes, qrule)
    fe = setup_test_mesh(6, shapes, qrule, 0.0, 1.0)

    # Set up globals and assemble system
    R = zeros(fe.nactive)
    K = zeros(fe.nactive, fe.nactive)
    assemble!(fe, R, K)

    # Factor, solve, and update
    update_U!(fe, K\R)

    # Check linear interpolation
    @test fe.U ≈ fe.mesh.X
end

test_fem1(Shapes1dP1(), GaussRule1dv(2))
test_fem1(Shapes1dP2(), GaussRule1dv(3))
test_fem1(Shapes1dP3(), GaussRule1dv(4))
