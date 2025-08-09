using LinearAlgebra
using Printf

#ldoc on
##
# # Mesh geometry
#
# A mesh consists of an array of nodes locations $x_j \in
# \mathbb{R}^d$ and an element connectivity array with `elt[i,j]`
# giving the node number for the $i$th node of the $j$th element.
#
# Each element represents a subset of $\Omega_e \subset \mathbb{R}^d$
# that is the image of a reference domain $\Omega_0 \subset
# \mathbb{R}^d$ under a mapping
# $$
#   \chi(\xi) = \sum_{i=1}^{m} N^e_i(\xi) x_i
# $$
# where $x_1, \ldots, x_m$ are the $m$ element node positions.  The
# functions $N^e_i$ are nodal basis functions (or Lagrange basis
# functions, or cardinal functions) for an interpolation set $\xi_1,
# \ldots, \xi_m \in \Omega_0$; that is $N_i(\xi_j) = \delta_{ij}$.
# The reference domain nodes $\xi_i$ are typically placed at corners
# or on edges of the reference domain, and their images are at
# corresponding locations in $\Omega_e$.
#
# When the same set of nodal basis functions (also called nodal shape
# functions in a finite element setting) are used both for defining
# the geometry and for approximating a PDE solution on $\Omega$, we
# call this method of describing the geometry an *isoparametric* map.
#
# We generally want our mappings describing the geometry to be
# *positively oriented*: that is, the map $\chi$ should be invertible
# and have positive Jacobian determinant over all of $\Omega_0$.
# This puts some restrictions on the spatial positions of the nodes;
# for example, if the interpolation nodes appear in counterclockwise
# order in the reference domain $\Omega_0$, then the corresponding
# spatial nodes in $\Omega_e$ should also appear in counterclockwise.
# order.
    
struct Mesh{T}
    shapes :: T             # Shape function interface
    X   :: Matrix{Float64}  # Node positions
    elt :: Matrix{Int}      # Connectivity
end

Mesh(shapes, numnp :: Integer, numelt :: Integer) =
    Mesh(shapes, zeros(dshapes(shapes), numnp),
         zeros(Int, nshapes(shapes), numelt))

##
# ## Block meshers
#
# One *can* allocate objects and then work out the node positions and
# element connectivity by hand (or with an external program).  But in
# many cases, a simpler option is to programatically generate a mesh
# that covers a simple domain (e.g. a block) and then map the
# locations of the nodes.  One can construct more complex meshes by
# combining this with a "tie" operation that merges the identity of
# nodes in the same location, but we will not bother with tied meshes
# for now.
#
# The simplest mesher creates a 1D mesh on an interval $[a,b]$.
# We allow elements of order 1-3.

function mesh_create1d(numelt, shapes, a=0.0, b=1.0)
    nen   = nshapes(shapes)
    numnp = numelt * (nen-1) + 1
    mesh = Mesh(shapes, numnp, numelt)
    mesh.X[:] .= range(a, b, length=numnp)
    for j = 1:numelt
        mesh.elt[:,j] .= (j-1)*(nen-1) .+ (1:nen)
    end
    mesh
end

##
# Things are more complicated in 2D, and we have distinct mesh
# generation routines for the different types of shape functions
# described in the `shapes` module.  Each of these generates a mesh
# of the region $[0,1]^2$ with `nex`-by-`ney` elements.

function mesh_block2d_P1(nex, ney)
    nx, ny = nex+1, ney+1
    mesh = Mesh(Shapes2dP1(), nx*ny, nex*ney)

    # Set up nodes (row-by-row, SW to NE)
    for iy = 1:ny
        for ix = 1:nx
            i = ix + (iy-1)*nx
            mesh.X[:,i] .= ( (ix-1)/(nx-1), (iy-1)/(ny-1) )
        end
    end

    # Set up element connectivity
    for iy = 1:ney
        for ix=1:nex
            i = ix + (iy-1)*nex
            i_sw = ix + (iy-1)*(nex+1)
            mesh.elt[:,i] .= (i_sw,
                              i_sw + 1,
                              i_sw + 1 + nex+1,
                              i_sw + nex+1)
        end
    end

    mesh
end

function mesh_block2d_P2(nex, ney)
    nx, ny = 2*nex+1, 2*ney+1
    mesh = Mesh(Shapes2dP2(), nx*ny, nex*ney)

    # Set up nodes (row-by-row, SW to NE)    
    for iy = 1:ny
        for ix = 1:nx
            i = ix + (iy-1)*nx
            mesh.X[:,i] .= ( (ix-1)/(nx-1), (iy-1)/(ny-1) )
        end
    end

    # Set up element connectivity
    for iy = 1:ney
        for ix = 1:nex
            i = ix + (iy-1)*nex
            i_sw = 2*(ix-1) + 2*(iy-1)*nx + 1
            mesh.elt[:,i] .= (i_sw,
                              i_sw + 1,
                              i_sw + 2,
                              i_sw + 2 + nx,
                              i_sw + 2 + 2*nx,
                              i_sw + 1 + 2*nx,
                              i_sw +   + 2*nx,
                              i_sw +   + nx,
                              i_sw + 1 + nx)
        end
    end
    
    mesh
end

function mesh_block2d_S2(nex, ney)
    nx0, nx1 = 2*nex+1, nex+1  # Even/odd row sizes
    numnp = (ney+1)*nx0 + ney*nx1
    mesh = Mesh(Shapes2dS2(), numnp, nex*ney)

    # Set up nodes (row-by-row, SW to NE)
    for iy = 1:ney
        start = (iy-1)*(nx0+nx1)

        # Fill bottom row
        for ix = 1:nx0
            mesh.X[:,start+ix] .= ( (ix-1)/(nx0-1), (iy-1)/ney )
        end

        # Fill middle row
        start += nx0
        for ix = 1:nx1
            mesh.X[:,start+ix] .= ( (ix-1)/(nx1-1), (iy-0.5)/ney )
        end
    end
    
    # Fill top row
    start = ney*(nx0+nx1)
    for ix = 1:nx0
        mesh.X[:,start+ix] .= ( (ix-1)/(nx0-1), 1.0 )
    end

    # Set up element connectivity
    for iy = 1:ney
        for ix = 1:nex
            i = ix + (iy-1)*nex
            i_sw = 2*(ix-1) + (iy-1)*(nx0+nx1) + 1
            i_ww =   (ix-1) + (iy-1)*(nx0+nx1) + nx0 + 1
            i_nw = 2*(ix-1) + (iy-1)*(nx0+nx1) + nx0 + nx1 + 1
            mesh.elt[:,i] .= (i_sw,
                              i_sw + 1,
                              i_sw + 2,
                              i_ww + 1,
                              i_nw + 2,
                              i_nw + 1,
                              i_nw,
                              i_ww)
        end
    end

    mesh
end

function mesh_block2d_T1(nex, ney)
    nx, ny = nex+1, ney+1
    mesh = Mesh(Shapes2dT1(), nx*ny, 2*nex*ney)

    # Set up nodes (row-by-row, SW to NE)
    for iy = 1:ney+1
        for ix = 1:nex+1
            i = ix + (iy-1)*(nex+1)
            mesh.X[:,i] .= ( (ix-1)/(nx-1), (iy-1)/(ny-1) )
        end
    end

    # Set up element connectivity
    for iy = 1:ney
        for ix = 1:nex
            i = ix + (iy-1)*nex;
            i_sw = ix + (iy-1)*(nex+1);

            # Two triangles makes a square
            mesh.elt[:,2*i-1] .= (i_sw, i_sw + 1, i_sw + nex+1)
            mesh.elt[:,2*i  ] .= (i_sw + nex+1, i_sw + 1, i_sw + 1 + nex+1)
        end
    end

    mesh
end

##
# ## Reference-to-spatial conversions
#
# Given a mesh and a point in a reference geometry (given by an
# element identifier `eltid` and coordinates `xref` in the element's
# reference domain), we would like to be able to compute spatial
# quantities (the shape functions, their spatial derivatives, and the
# Jacobian of the reference to spatial map).  The Jacobian matrix
# is in LU-factored form.
#
# NB: This gets used a lot -- at every mesh point and at every stage.
# So we probably don't want to do very much dynamic allocation.  We
# are not so careful about this now; this should perhaps be revisited.

function mesh_to_spatial!(mesh, eltid, x, J)
    s = mesh.shapes(x)
    e = view(mesh.elt,:,eltid)
    mul!(x, view(mesh.X,:,e), s.N, 1.0, 0.0)
    mul!(J, view(mesh.X,:,e), s.dN, 1.0, 0.0)
    F = lu(J)  # NB: May want to do this in place
    rdiv!(s.dN, F)
    x, F, det(F)
end

function mesh_to_spatial(mesh, eltid, xref)
    d = dshapes(mesh.shapes)
    J = zeros(d,d)
    mesh_to_spatial!(mesh, eltid, copy(xref), J)
end

##
# ## Mesh output
# 
# For debugging, it is helpful to be able to print out all or part of
# the mesh geometry.  We mostly care about this for looking at small
# meshes.

function mesh_print_nodes(mesh)
    @printf("\nNodal Positions:\n")
    @printf("   ID ")
    for j = 1:dshapes(mesh.shapes)
        @printf("     X%d", j)
    end
    @printf("\n")
    for i = 1:size(mesh.X,2)
        @printf("%3d : ", i)
        for j = 1:dshapes(mesh.shapes)
            @printf(" %6.2g", mesh.X[j,i])
        end
        @printf("\n")
    end
end

function mesh_print_elt(mesh)
    @printf("\nElement connectivity:\n")
    for i = 1:size(mesh.elt,2)
        @printf("% 3d :", i)
        for j = 1:size(mesh.elt,1)
            @printf("  % 3d", mesh.elt[j,i])
        end
        @printf("\n")
    end
end

function mesh_print(mesh)
    mesh_print_nodes(mesh)
    mesh_print_elt(mesh)
end
