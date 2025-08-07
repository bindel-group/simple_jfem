mutable struct FEMProblem{T}

    # Mesh data
    mesh :: Mesh

    # Element type (NB: can generalize with multiple types)
    etype :: T

    # Storage for fields
    U  :: Matrix{Float64}  # Global soln values (ndof-by-numnp)
    F  :: Matrix{Float64}  # Global force values (ndof-by-numnp)
    id :: Matrix{Integer}  # Global to reduced ID map (ndof-by-nump)

    # Dimensions
    ndof :: Integer
    nactive :: Integer

end

function FEMProblem(mesh, etype, ndof)
    numnp = size(mesh.X,2)
    nactive = numnp * ndof
    U = zeros(ndof, numnp)
    F = zeros(ndof, numnp)
    id = zeros(Integer, ndof, numnp)
    FEMProblem(mesh, etype, U, F, id, ndof, nactive)
end

function assign_ids(fe :: FEMProblem)
    nactive = 0
    for j = 1:size(fe.mesh.X,2)
        for i = 1:fe.ndof
            if id[i,j] >= 0
                nactive += 1
                id[i] = nactive
            end
        end
    end
    fe.nactive = nactive
end

function update_U(fe :: FEMProblem, du_red)
    for j = 1:size(fe.mesh.X,2)
        for i = 1:fe.ndof
            if id[i,j] > 0
                U[i,j] -= du_red[id[i,j]]
            end
        end
    end
end

function set_load(fe :: FEMProblem, f :: Function)
    for i = 1:size(fe.mesh.X,2)
        f(view(fe.mesh.X,:,i), view(fe.F,:,i))
    end
end

function assemble(fe :: FEMProblem, R :: Vector, K)
    nlocal = nshapes(mesh.shapes) * fe.ndof
    Re = zeros(nlocal)
    Ke = zeros(nlocal,nlocal)
    ids = zeros(Integer, nlocal)
    id[elt]
    for i = 1:size(fe.mesh.elt,2)
        ids[:] .= view(fe.id,:,view(elt,:,j))
        element_dR(etype, fe, i, Re, Ke)
        assemble(R, Re, ids)
        assemble(K, Ke, ids)
    end
    R, K
end

function assemble(fe :: FEMProblem, R :: Vector)
    nlocal = nshapes(mesh.shapes) * fe.ndof
    Re = zeros(nlocal)
    ids = zeros(Integer, nlocal)
    id[elt]
    for i = 1:size(fe.mesh.elt,2)
        ids[:] .= view(fe.id,:,view(elt,:,j))
        element_dR(etype, fe, i, Re)
        assemble(R, Re, ids)
    end
    R
end

function assemble(fe :: FEMProblem, K)
    nlocal = nshapes(mesh.shapes) * fe.ndof
    Ke = zeros(nlocal,nlocal)
    ids = zeros(Integer, nlocal)
    id[elt]
    for i = 1:size(fe.mesh.elt,2)
        ids[:] .= view(fe.id,:,view(elt,:,j))
        element_dR(etype, fe, i, Ke)
        assemble(K, Ke, ids)
    end
    K
end

function fem_print(fe :: FEMProblem)
    @printf("\nNodal information:\n")
    @printf("       ID ")
    for j = 1:dshapes(fe.mesh.shapes)
        @printf("     X%d", j)
    end
    for j = 1:fe.ndof
        @printf("     U%d", j)
    end
    for j = 1:fe.ndof
        @printf("     F%d", j)
    end
    @printf("\n")
    for i = 1:size(fe.mesh.X,2)
        @printf("%3d : ", i)
        for j = 1:ndof
            @printf("% 3d ", fe.id[i])
        end
        for j = 1:dshapes(fe.mesh.shapes)
            @printf(" %6.2g", fe.mesh.X[j,i])
        end
        for j = 1:fe.ndof
            @printf(" % 6.2g", fe.U[j,i])
        end
        for j = 1:fe.ndof
            @printf(" % 6.2g", fe.F[j,i])
        end
        @printf("\n");
    end
    mesh_print_elt(fe.mesh)
end
