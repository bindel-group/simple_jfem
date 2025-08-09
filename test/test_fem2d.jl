include("../src/quadrules.jl")
include("../src/shapes.jl")
include("../src/mesh.jl")
include("../src/assemble.jl")
include("../src/fem.jl")
include("../src/element.jl")

using Test

@testset "2D FE test" begin

# Set up mesh on [0,1]^2 with Dirichlet BC
mesh = mesh_block2d_P1(2, 2)
fe = FEMProblem(mesh, PoissonElt(), GaussRule2d(3), 1)

# Move midpoint to off center (patch test)
mesh.X[:,5] .= (0.6; 0.4)

# BC at x = 0 and 1
function dirichlet_bc(x, ids, u)
    if x[1] == 0.0
        ids[1] = -1
        u[1] = 0.0
    elseif x[1] == 1.0
        ids[1] = -1
        u[1] = 1.0
    end
end

set_dirichlet!(fe, dirichlet_bc)

# Set up and solve
assign_ids!(fe)
R = zeros(fe.nactive)
K = zeros(fe.nactive, fe.nactive)
assemble!(fe, R, K)
update_U!(fe, K\R)

# Check vs reference solution
@test fe.U[1,:] ≈ fe.mesh.X[1,:]

end
