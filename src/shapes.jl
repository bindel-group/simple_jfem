function shapes1dP1(x)
    0.5*(1-x),
    0.5*(1+x)
end

function dshapes1dP1(x)
    -0.5,
    0.5
end

function shapes1dP1!(N, dN, xx)
    x = xx[1]
    N[:]  .= shapes1dP1(x)
    dN[:] .= dshapes1dP1(x)
end

function shapes1dP2(x)
    -0.5*(1-x)*x,
         (1-x)*  (1+x),
     0.5*      x*(1+x)
end

function dshapes1dP2(x)
    -0.5*(1-2*x),
    -2*x,
     0.5*(1+2*x)
end

function shapes1dP2!(N, dN, xx)
    x = xx[1]
    N[:]  .= shapes1dP2(x)
    dN[:] .= dshapes1dP2(x)
end

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

function shapes1dP3!(N, dN, xx)
    x = xx[1]
    N[:]  .= shapes1dP3(x)
    dN[:] .= dshapes1dP3(x)
end

function shapes2dP1!(N, dN, xx)
    Nx1,  Nx2  = shapes1dP1(xx[1])
    Ny1,  Ny2  = shapes1dP1(xx[2])
    dNx1, dNx2 = dshapes1dP1(xx[1])
    dNy1, dNy2 = dshapes1dP1(xx[2])
    N[:]  .= ( Nx1* Ny1,  Nx2* Ny1,  Nx2* Ny2,  Nx1* Ny2)
    dN[:] .= (dNx1* Ny1, dNx2* Ny1, dNx2* Ny2, dNx1* Ny2,
               Nx1*dNy1,  Nx2*dNy1,  Nx2*dNy2,  Nx1*dNy2)              
end

function shapes2dP2!(N, dN, xx)
    Nx1,  Nx2,  Nx3  = shapes1dP2(xx[1])
    Ny1,  Ny2,  Ny3  = shapes1dP2(xx[2])
    dNx1, dNx2, dNx3 = dshapes1dP2(xx[1])
    dNy1, dNy2, dNy3 = dshapes1dP2(xx[2])
    N[:] .=  ( Nx1* Ny1,  Nx2* Ny1,  Nx3* Ny1,
               Nx3* Ny2,  Nx3* Ny3,  Nx2* Ny3,
               Nx1* Ny3,  Nx1* Ny2,  Nx2* Ny2)
    dN[:] .= (dNx1* Ny1, dNx2* Ny1, dNx3* Ny1,
              dNx3* Ny2, dNx3* Ny3, dNx2* Ny3,
              dNx1* Ny3, dNx1* Ny2, dNx2* Ny2,
               Nx1*dNy1,  Nx2*dNy1,  Nx3*dNy1,
               Nx3*dNy2,  Nx3*dNy3,  Nx2*dNy3,
               Nx1*dNy3,  Nx1*dNy2,  Nx2*dNy2)
end

function shapes2dS2!(N, dN, xx)
    Nx1,  Nx2,  Nx3  = shapes1dP2(xx[1])
    Ny1,  Ny2,  Ny3  = shapes1dP2(xx[2])
    dNx1, dNx2, dNx3 = dshapes1dP2(xx[1])
    dNy1, dNy2, dNy3 = dshapes1dP2(xx[2])
    N[:] .=  ( Nx1* Ny1,  Nx2* Ny1,  Nx3* Ny1,
               Nx3* Ny2,  Nx3* Ny3,  Nx2* Ny3,
               Nx1* Ny3,  Nx1* Ny2)
    dN[:] .= (dNx1* Ny1, dNx2* Ny1, dNx3* Ny1,
              dNx3* Ny2, dNx3* Ny3, dNx2* Ny3,
              dNx1* Ny3, dNx1* Ny2,
               Nx1*dNy1,  Nx2*dNy1,  Nx3*dNy1,
               Nx3*dNy2,  Nx3*dNy3,  Nx2*dNy3,
               Nx1*dNy3,  Nx1*dNy2)
end

function shapes2dT1!(N, dN, xx)
    N[:]  .= (1.0-xx[1]-xx[2], xx[1], xx[2])
    dN[:] .= (-1.0, 1.0, 0.0,
              -1.0, 0.0, 1.0)
end
