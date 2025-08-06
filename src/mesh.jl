using LinearAlgebra
using Printf

struct Mesh{T}
    shapes :: T
    X   :: Matrix{Float64}
    elt :: Matrix{Integer}
end

Mesh(shapes, numnp :: Integer, numelt :: Integer) =
    Mesh(shapes, zeros(dshapes(shapes), numnp),
         zeros(Integer, nshapes(shapes), numelt))

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

function mesh_block2d_P1(nex, ney)
    nx, ny = nex+1, ney+1
    mesh = Mesh(Shapes2dP1(), nx*ny, nex*ney)

    # Set up nodes (row-by-row, SW to NE)
    for iy = 1:nx
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
            mesh.elt[:,i] = (i_sw,
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
            mesh[:,i] .= ( (ix-1)/(nx-1), (iy-1)/(ny-1) )
        end
    end

    # Set up element connectivity
    for iy = 1:ney
        for ix = 1:nex
            i = ix + (iy-1)*nex
            i_sw = 2*(ix-1) + 2*(iy-1)*nx + 1
            mesh.elt[:,i] = (i_sw,
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
    mesh = Mesh(Shape2dS2(), numnp, nex*ney)

    # Set up nodes (row-by-row, SW to NE)
    for iy = 1:ney
        start = (iy-1)*(nx0+nx1)

        # Fill bottom row
        for ix = 1:nx0
            mesh.X[:,start+ix] = ( (ix-1)/(nx0-1), (iy-1)/ney )
        end

        # Fill middle row
        start += nx0
        for ix = 1:nx1
            mesh.X[:,start+ix] = ( (ix-1)/(nx1-1), (iy-0.5)/ney )
        end
    end

    # Set up element connectivity
    for iy = 1:ney
        for ix = 1:nex
            i = ix + (iy-1)*nex
            i_sw = 2*(ix-1) + (iy-1)*(nx0+nx1) + 1
            i_ww =   (ix-1) + (iy-1)*(nx0+nx1) + nx0 + 1
            i_nw = 2*(ix-1) + (iy-1)*(nx0+nx1) + nx0 + nx1 + 1
            mesh.elt[:,i] = (i_sw,
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
            mesh.X[:,i] = ( (ix-1)/(nx-1), (iy-1)/(ny-1) )
        end
    end

    # Set up element connectivity
    for iy = 1:ney
        for ix = 1:nex
            i = ix + (iy-1)*nex;
            i_sw = ix + (iy-1)*(nex+1);

            # Two triangles makes a square
            mesh.elt[:,2*i-1] = (i_sw, i_sw + 1, i_sw + nex+1)
            mesh.elt[:,2*i  ] = (i_sw + nex+1, i_sw + 1, i_sw + 1 + nex+1)
        end
    end

    mesh
end

function mesh_to_spatial!(mesh, eltid, xref, F)
    s = mesh.shapes(xref)
    xref[:] .= 0.0
    F.factors[:] .= 0.0
    for k = 1:nshapes(s)
        j = mesh.elt[k,eltid]
        Nj = N[j]
        xj = @view mesh.X[:,j]
        dNj = @view dN[j,:]
        @. xref[:] += xj * s.N[j]
        mul!(F.factors, xj, dNj, 1.0, 1.0)
    end
    LAPACK.getrf!(F.factors, F.ipiv)
    rdiv!(s.dN, F)
    x, F, det(F)
end

function mesh_to_spatial(mesh, eltid, xref)
    F = LU(zeros(d,d), zeros(Integer,d), zero(LAPACK.BlasInt))
    mesh_to_spatial!(mesh, eltid, copy(xref), F)
end

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
            @printf(" %6.2g", X[i,j])
        end
        @printf("\n")
    end
end

function mesh_print_elt(mesh)
    @printf("\nElement connectivity:\n")
    for i = 1:size(mesh.elts,2)
        @printf("% 3d :", i)
        for j = 1:size(mesh.elts,1)
            @printf("  % 3d", mesh.elt[j,i])
        end
        @printf("\n")
    end
end

function mesh_print(mesh)
    mesh_print_nodes(mesh)
    mesh_print_elt(mesh)
end
