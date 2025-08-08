struct Poisson1dElt end
struct Poisson2dElt end

# Dispatch to specialized version based on element type
element_dR!(fe :: FEMProblem, eltid, Re, Ke) =
    element_dR!(fe.etype, fe, eltid, Re, Ke)

# The Poisson 1D and 2D codes are essentially identical
element_dR!(:: Poisson1dElt, fe :: FEMProblem, eltid, Re, Ke) =
    poisson_element_dR(fe, eltid, Re, Ke)
element_dR!(:: Poisson2dElt, fe :: FEMProblem, eltid, Re, Ke) =
    poisson_element_dR(fe, eltid, Re, Ke)

function poisson_element_dR(fe :: FEMProblem, eltid, Re, Ke)
    s = fe.mesh.shapes
    eltj = view(fe.mesh.elt,:,eltid)
    for (xi,wt) in fe.qrule
        x, f, detJ = mesh_to_spatial(fe.mesh, eltid, xi)
        wt *= detJ
        du = (view(fe.U,1,eltj)' * s.dN)'
        fx =  view(fe.F,1,eltj)' * s.N
        Re[:] += (s.dN*du - s.N*fx) * wt
        mul!(Ke, s.dN, s.dN', wt, 1.0)
    end
end
