#ldoc on
##
# # Finite element mesh
#
# My finite element mesh data structure is informed by lots of
# old Fortran codes, and mostly is a big pile of arrays.
# Specifically, we have the nodal arrays:
#
# - `U`: Global array of solution values, *including* those that are
#   determined by Dirichlet boundary conditions.  Column $j$
#   represents the unknowns at node $j$ in the mesh.
# - `F`: Global array of load values (right hand side evaluations of
#   the forcing function in Poisson, for example; but Neumann boundary
#   conditions can also contribute to `F`).
# - `id`: Indices of solution values in a reduced solution vector.
#   One column per node, with the same dimensions as `U` (and `F`), so
#   that `ureduced[id[i,j]]` corresponds to `U[i,j]` when `id[i,j]` is
#   nonnegative.  The reduced solution vector contains only those
#   variables that are not constrained a priori by boundary
#   conditions; we mark the latter with negative entries in the `id`
#   array.
# 
# In addition, we keep a mesh, an element type, and a quadrature rule.
# Note that for the moment, we are assuming only one element type per
# problem; we could have a separate array of element types (one per
# element) if we wanted more flexibility.

mutable struct FEMProblem{T,S}

    # Mesh data
    mesh :: Mesh

    # Element type (NB: can generalize with multiple types)
    etype :: T

    # Quadrature rule
    qrule :: S

    # Storage for fields
    U  :: Matrix{Float64}  # Global soln values (ndof-by-numnp)
    F  :: Matrix{Float64}  # Global force values (ndof-by-numnp)
    id :: Matrix{Integer}  # Global to reduced ID map (ndof-by-nump)

    # Dimensions
    ndof :: Integer
    nactive :: Integer

end

function FEMProblem(mesh, etype, qrule, ndof)
    numnp = size(mesh.X,2)
    nactive = numnp * ndof
    U = zeros(ndof, numnp)
    F = zeros(ndof, numnp)
    id = zeros(Integer, ndof, numnp)
    FEMProblem(mesh, etype, qrule, U, F, id, ndof, nactive)
end

##
# ## Index setup
#
# The `assign_ids!` function sets up the `id` array.  On input, the
# `id` entries should be initialized so that boundary values are
# marked with negative numbers, and everything else is non-negative.
# On output, entries of `id` for variables not subject to essential
# boundary conditions will be assigned indices from 1 to `nactive`
# (and `nactive` will be updated appropriately).

function assign_ids!(fe :: FEMProblem)
    nactive = 0
    for j = 1:size(fe.mesh.X,2)
        for i = 1:fe.ndof
            if fe.id[i,j] >= 0
                nactive += 1
                fe.id[i,j] = nactive
            end
        end
    end
    fe.nactive = nactive
end

##
# ## Solution updates
#
# The `update_U!` function applies an update to the internal state.
# That we compute `U[i,j] -= du_red[id[i,j]]` for `id[i,j] > 0`.  If
# the update comes from $K^{-1} R$ where $K$ is the reduced tangent
# and $R$ the reduced residual, then applying the update will exactly
# solve the equation in the linear PDE case.  However, we can also
# appy approximate updates (e.g. with an inexact solver for $K$), and
# the same framework works for Newton iterations for nonlinear
# problems.

function update_U!(fe :: FEMProblem, du_red)
    for j = 1:size(fe.mesh.X,2)
        for i = 1:fe.ndof
            if fe.id[i,j] > 0
                fe.U[i,j] -= du_red[fe.id[i,j]]
            end
        end
    end
end

##
# ## Loads and BCs
#
# The `set_load!` function and `set_dirichlet!` function update the
# forcing array `F` and the boundary data entries and `id` array markers
# in `U` and `id`, respectively.

function set_load!(fe :: FEMProblem, f :: Function)
    for i = 1:size(fe.mesh.X,2)
        f(view(fe.mesh.X,:,i), view(fe.F,:,i))
    end
end

function set_dirichlet!(fe :: FEMProblem, f :: Function)
    for i = 1:size(fe.mesh.X,2)
        f(view(fe.mesh.X,:,i), view(fe.id,:,i), view(fe.U,:,i))
    end
end

##
# ## Assembly
#
# The assembly loops iterate through the elements and produce a global
# residual and tangent stiffness based on the current solution state.

function assemble!(fe :: FEMProblem, R :: Vector, K)
    nlocal = nshapes(fe.mesh.shapes) * fe.ndof
    Re = zeros(nlocal)
    Ke = zeros(nlocal,nlocal)
    ids = zeros(Integer, nlocal)
    clear!(R)
    clear!(K)
    for i = 1:size(fe.mesh.elt,2)
        ids[:] .= reshape(view(fe.id,:,view(fe.mesh.elt,:,i)),:)
        Re[:] .= 0.0
        Ke[:] .= 0.0
        element_dR!(fe, i, Re, Ke)
        assemble_add!(R, Re, ids)
        assemble_add!(K, Ke, ids)
    end
    R, K
end

##
# ## Debugging printer

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
