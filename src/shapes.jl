#ldoc on
##
# # Shape functions
#
# A *shape function* on a reference domain is a basis function used
# for interpolation on that domain.  We will generally use Lagrange
# shape functions (also called nodal shape functions), which are one
# at one nodal point in a reference domain and zero at the others.
# We want to be able to compute both the values of all shape functions
# at a point in the domain and also their derivatives (stored as
# a matrix with $d$ columns for a $d$-dimensional reference domain).
#
# ## 1D building blocks
#
# A useful building block is the Lagrange polynomials through 2, 3, or
# 4 equispaced points, which span the linear (P1), quadratic (P2), and
# cubic (P3) polynomial spaces on $[-1, 1]$, respectively.

shapes1dP1(x) = (0.5*(1-x), 0.5*(1+x))
dshapes1dP1(x) =(-0.5,      0.5)

shapes1dP2(x) = (-0.5*(1-x)*x, (1-x)*(1+x), 0.5*x*(1+x))
dshapes1dP2(x) =(-0.5*(1-2*x), -2*x,        0.5*(1+2*x))

function shapes1dP3(x)
    (-1.0/16) * (1-x)*(1-3*x)*(1+3*x),
    ( 9.0/16) * (1-x)*(1-3*x)*(1+x),
    ( 9.0/16) * (1-x)*(1+3*x)*(1+x),
    (-1.0/16) * (1-3*x)*(1+3*x)*(1+x)
end

function dshapes1dP3(x)
    1.0/16 * ( 1+x*( 18+x*-27)),
    9.0/16 * (-3+x*(-2+x* 9)),
    9.0/16 * ( 3+x*(-2+x*-9)),
    1.0/16 * (-1+x*( 18+x* 27))
end

##
# ## Shape class macro
# 
# We make a general structure type for accessing shape functions.  For
# a given set of shapes, we keep some storage for the shape functions
# (`N`) and their derivatives (`dN`).  Because we are defining the same
# types of methods for each instance, we use Julia's macro facility to
# help write the routines for us.

macro make_shape(classname, nfun, dim, body)
    quote
        struct $(esc(classname))
            N  :: Vector{Float64}  # Shape functions (nshapes)
            dN :: Matrix{Float64}  # Derivatives (nshapes-by-dshapes)
        end

        # Default constructor
        $(esc(classname))() =
            $(esc(classname))(zeros($nfun), zeros($nfun, $dim))

        # Get the number of shapes and dimension of the space
        $(esc(:nshapes))( :: $(esc(classname))) = $nfun
        $(esc(:dshapes))( :: $(esc(classname))) = $dim

        # Routine for actual computation (looks like a call to the object)
        function (s :: $(esc(classname)))(xx)
            $body
            s
        end
    end
end

##
# ## 1D shape functions
#
# The 1D shape function classes just call the Lagrange polynomial
# support routines defined earlier.  For each, shape, we also define
# a `refnodes` array with the list of reference-domain locations for
# the nodes.

@make_shape Shapes1dP1 2 1 begin
    x = xx[1]
    s.N[:]  .= shapes1dP1(x)
    s.dN[:] .= dshapes1dP1(x)
end
refnodes(:: Shapes1dP1) = [-1.0 1.0]

@make_shape Shapes1dP2 3 1 begin
    x = xx[1]
    s.N[:]  .= shapes1dP2(x)
    s.dN[:] .= dshapes1dP2(x)
end
refnodes(:: Shapes1dP2) = [-1.0 0.0 1.0]

@make_shape Shapes1dP3 4 1 begin
    x = xx[1]
    s.N[:]  .= shapes1dP3(x)
    s.dN[:] .= dshapes1dP3(x)
end
refnodes(:: Shapes1dP3) = [-1.0 -1.0/3 1.0/3 1.0]

##
# ## 2D shape functions
#
# The P1 and P2 shapes are just tensor products of the 1D P1 and P2
# shapes.  The serendipity element S2 is like the P2, but doesn't
# keep the "bubble mode" shape function associated with the center node.

@make_shape Shapes2dP1 4 2 begin
    Nx1,  Nx2  = shapes1dP1(xx[1])
    Ny1,  Ny2  = shapes1dP1(xx[2])
    dNx1, dNx2 = dshapes1dP1(xx[1])
    dNy1, dNy2 = dshapes1dP1(xx[2])
    s.N[:]  .= ( Nx1* Ny1,  Nx2* Ny1,  Nx2* Ny2,  Nx1* Ny2)
    s.dN[:] .= (dNx1* Ny1, dNx2* Ny1, dNx2* Ny2, dNx1* Ny2,
                 Nx1*dNy1,  Nx2*dNy1,  Nx2*dNy2,  Nx1*dNy2)
end
refnodes(:: Shapes2dP1) =
    [-1.0   1.0   1.0  -1.0 ;
     -1.0  -1.0   1.0   1.0 ]

@make_shape Shapes2dP2 9 2 begin
    Nx1,  Nx2,  Nx3  = shapes1dP2(xx[1])
    Ny1,  Ny2,  Ny3  = shapes1dP2(xx[2])
    dNx1, dNx2, dNx3 = dshapes1dP2(xx[1])
    dNy1, dNy2, dNy3 = dshapes1dP2(xx[2])
    s.N[:] .=  ( Nx1* Ny1,  Nx2* Ny1,  Nx3* Ny1,
                 Nx3* Ny2,  Nx3* Ny3,  Nx2* Ny3,
                 Nx1* Ny3,  Nx1* Ny2,  Nx2* Ny2)
    s.dN[:] .= (dNx1* Ny1, dNx2* Ny1, dNx3* Ny1,
                dNx3* Ny2, dNx3* Ny3, dNx2* Ny3,
                dNx1* Ny3, dNx1* Ny2, dNx2* Ny2,
                 Nx1*dNy1,  Nx2*dNy1,  Nx3*dNy1,
                 Nx3*dNy2,  Nx3*dNy3,  Nx2*dNy3,
                 Nx1*dNy3,  Nx1*dNy2,  Nx2*dNy2)
end
refnodes(:: Shapes2dP2) =
    [-1.0   0.0   1.0   1.0   1.0   0.0  -1.0  -1.0   0.0 ;
     -1.0  -1.0  -1.0   0.0   1.0   1.0   1.0   0.0   0.0 ]

@make_shape Shapes2dS2 8 2 begin
    x, y = xx
    Nx1,  Nx2,  Nx3  = shapes1dP2(x)
    Ny1,  Ny2,  Ny3  = shapes1dP2(y)
    dNx1, dNx2, dNx3 = dshapes1dP2(x)
    dNy1, dNy2, dNy3 = dshapes1dP2(y)
    s.N[:] .=  ( Nx1* Ny1,  Nx2* Ny1,  Nx3* Ny1,
                 Nx3* Ny2,  Nx3* Ny3,  Nx2* Ny3,
                 Nx1* Ny3,  Nx1* Ny2)
    s.dN[:] .= (dNx1* Ny1, dNx2* Ny1, dNx3* Ny1,
                dNx3* Ny2, dNx3* Ny3, dNx2* Ny3,
                dNx1* Ny3, dNx1* Ny2,
                 Nx1*dNy1,  Nx2*dNy1,  Nx3*dNy1,
                 Nx3*dNy2,  Nx3*dNy3,  Nx2*dNy3,
                 Nx1*dNy3,  Nx1*dNy2)
end
refnodes(:: Shapes2dS2) =
    [-1.0   0.0   1.0   1.0   1.0   0.0  -1.0  -1.0   0.0 ;
     -1.0  -1.0  -1.0   0.0   1.0   1.0   1.0   0.0   0.0 ]

##
# ## Triangle shapes
#
# For now, we only have linear shape functions on a triangle.

@make_shape Shapes2dT1 3 2 begin
    x, y = xx
    s.N[:]  .= (1.0-x-y, x, y)
    s.dN[:] .= (-1.0, 1.0, 0.0,
                -1.0, 0.0, 1.0)
end
refnodes(:: Shapes2dT1) =
    [ 0.0 1.0 0.0 ;
      0.0 0.0 1.0 ]
