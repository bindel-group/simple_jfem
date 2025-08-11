#ldoc on
##
# # Elements
#
# Abstractly, for steady-state problems, we are finding 
# $u(x) = \sum_j N_j(x) u_j$ via an equation
# $$
#   R(u, N_i) = 0
# $$
# for all shape functions $N_i$ that are not associated with
# essential boundary conditions.  The element routines compute
# the contribution of one element to the residual $R$ and to 
# the tangent $\partial R/\partial u_j$.
# 
# Different types of equations demand different types of elements.
# Even for a single type of element, we may depend on things like
# PDE coefficients or choices of material parameters (as well as
# implementation details like the quadrature rule used for computing
# integrals).  An `element_t` object type keeps all this information
# together.  The `element_t` data type should be thought of as
# representing a *type* of element, and not one specific element;
# usually many elements share fundamentally the same data, differing
# only in which nodes they involve.  In the language of design patterns,
# this is an example of a "flyweight" pattern.
# 
# The main interface for an element is a method
# 
#     element_dR!(etype, fe, eltid, Re, Ke)
# 
# where `etype` is data for the element type, `fe` is a finite
# element mesh data structure, `eltid` is the index of the element in
# the mesh, and `Re` and `Ke` are pointers to storage for the element
# residual and tangent matrix contributions.  Either `Re` or `Ke` can
# be `nothing`, indicating that we don't need that output.
#
# Because we will usually use the element type associated with the
# finite element problem, we provide a convenience wrapper that fills
# in the `etype` argument to `element_dR!` with the problem `etype`.

element_dR!(fe :: FEMProblem, eltid, Re, Ke) =
    element_dR!(fe.etype, fe, eltid, Re, Ke)

##
# ## Poisson element
#
# Right now, we only have one element type, corresponding to a Poisson
# problem, written in weak form as
# $$
#  R(u, N_i) =
#  \int_{\Omega} \left(
#  \nabla N_i(x) \cdot \nabla u(x) -
#  N_i(x) f(x) \right) \, d\Omega(x).
# $$
# There are no PDE coefficients or other special parameters to keep
# track of for this element type.  This is also a simple enough
# problem that it's easy to write dimension-independent code -- no
# need for distinguishing 1D from 2D elements.
#
# Note: We are a little careful to avoid dynamic allocation inside
# `element_dR!`, as this is called many times.

struct PoissonElt
    du :: Vector{Float64}  # Workspace for element computations
end

PoissonElt(d :: Integer) = PoissonElt(zeros(d))

function element_dR!(e :: PoissonElt, fe :: FEMProblem, eltid, Re, Ke)
    s, J, ipiv = fe.mesh.shapes, fe.mesh.J, fe.mesh.ipiv
    eltj = view(fe.mesh.elt,:,eltid)
    X = view(fe.mesh.X,:,eltj)
    U = view(fe.U,1,eltj)
    F = view(fe.F,1,eltj)
    du = e.du
    for (x,wt) in IsoMappedRule(fe.qrule, s, X, true, J, ipiv)
        fx =  F' * s.N
        mul!(du, s.dN, U)
        mul!(Re, s.dN', du,   wt, 1.0)
        mul!(Re, s.N,   fx,  -wt, 1.0)
        mul!(Ke, s.dN', s.dN, wt, 1.0)
    end
end
