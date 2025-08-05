# -- 1D shapes and derivatives

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

# -- Shape function class setup

macro make_shape(classname, nfun, dim, body)
    quote
        struct $(esc(classname))
            N  :: Vector{Float64}
            dN :: Matrix{Float64}
        end

        $(esc(classname))() =
            $(esc(classname))(zeros($nfun), zeros($nfun, $dim))

        $(esc(:nshapes))( :: $(esc(classname))) = $nfun
        $(esc(:dshapes))( :: $(esc(classname))) = $dim

        function (s :: $(esc(classname)))(xx)
            $body
            s
        end
    end
end

# -- Shape function class definitions

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

@make_shape Shapes2dT1 3 2 begin
    x, y = xx
    s.N[:]  .= (1.0-x-y, x, y)
    s.dN[:] .= (-1.0, 1.0, 0.0,
                -1.0, 0.0, 1.0)
end
refnodes(:: Shapes2dT1) =
    [ 0.0 1.0 0.0 ;
      0.0 0.0 1.0 ]
