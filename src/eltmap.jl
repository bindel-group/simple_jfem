using LinearAlgebra

#ldoc on
##
# # Mapped subdomains
#
# Our shape functions and quadrature rules are defined on simple
# reference domains.  To accommodate more complicated domains, we
# consider mappings $\chi : \mathbb{R}^d \rightarrow \mathbb{R}^d$
# from a reference domain $\Omega_0$ into a spatial domain $\Omega_e$
# (with $\chi(\Omega_0) = \Omega_e$).  We will typically write
# $x = \chi(X)$ for generic points in the spatial and reference domains.
#
# ## Mapping functions and derivatives
# 
# With the mapping between domains in hand, We can identify functions
# $f^0$ in the reference domain with functions $f^e$ in the spatial
# domain:
# $$
#   f^e(x) = (f^0 \circ \chi^{-1})(x) = f^0(X)
# $$
# and the derivative mapping is given by
# $$
#   {f^e}'(x) = {f^0}'(X) J(X)^{-1}
# $$
# where $J(X) = \partial \chi/\partial X$ is the Jacobian of the
# reference-to-spatial mapping.  In terms of gradients, we can
# write this as
# $$
#   \nabla_x f^e(x) = J(X)^{-T} \, \nabla_X f^0(X).
# $$
#
# Because we expect to be computing these types of transformations
# a lot, we would rather not constantly allocate and deallocate space
# for LU factorization objects.  Therefore, we define helper functions
# that directly call the LAPACK interfaces for LU factorization of $J$
# and subsequent solves with $J^T$ (for a generic gradient $g$ or for
# the shape function gradients stored in a `ShapeFuns` struct).

r2s_factor!(J, ipiv) = LAPACK.getrf!(J, ipiv)
r2s_gradients!(J, ipiv, g :: AbstractMatrix) = LAPACK.getrs!('T', J, ipiv, g)
r2s_gradients!(J, ipiv, s :: ShapeFuns) = r2s_gradients!(J, ipiv, s.dN)

##
# We also define a function to get the determinant out of the
# LU factorization of $J$.

function r2s_det(J, ipiv)
    detJ = 1.0
    for k = 1:size(J,2)
        detJ *= J[k,k]
        if ipiv[k] != k
            detJ = -detJ
        end
    end
    detJ
end

##
# ## Mapped quadrature
#
# We can also express integrals over $\Omega_e$ in terms of integrals
# over $\Omega_0$ via the usual change of variables formula:
# $$
#   \int_{\Omega_e} f^e(x) \, dx =
#   \int_{\Omega_0} f^0(X) |\det J(X)| \, dX.
# $$
# We will generally assume that $\chi$ is *positively-oriented*, that is,
# $\det J$ is positive over all points within $\Omega_0$.  Under this
# assumption, we can convert quadrature rules over the reference domain
# to mapped quadrature rules:
# $$
#   \int_{\Omega_e} f^e(x) \, dx \approx
#   \sum_{j=1}^p f^e(\chi(\xi_j)) \det J(\xi_j) w_j.
# $$

struct MappedRule{T,M} <: QuadratureRule where {T <: QuadratureRule}
    base_rule :: T                # Reference domain rule
    chi!      :: M                # Mapping
    J         :: Matrix{Float64}  # Space for Jacobian/LU
    ipiv      :: Vector{Int}      # Pivots for Jacobian LU
end

MappedRule(base_rule, chi!) =
    MappedRule(base_rule, chi!,
               zeros(quad_dim(base_rule), quad_dim(base_rule)),
               zeros(Int, quad_dim(base_rule)))

r2s_gradients!(q :: MappedRule, g) = r2s_gradients!(q.J, q.ipiv, g)

quad_npoints(q :: MappedRule) = quad_npoints(q.base_rule)
quad_dim(q :: MappedRule) = quad_dim(q.base_rule)

function quad_point(q :: MappedRule, i)
    x = quad_point(q.base_rule, i)
    q.chi!(x, q.J)
    r2s_factor!(q.J, q.ipiv)
    x
end

##
# By default, `quad_point` is called before `quad_weight`, and so we
# do not need to re-evaluate the mapping and its Jacobian

function quad_weight(q :: MappedRule, i; map = false)
    if map quad_point(q, i) end
    quad_weight(q.base_rule, i) * r2s_det(q.J, q.ipiv)
end

##
# ## Iso-parametric mapping
#
# In some cases, the mapping function $\chi$ may be hand-crafted.
# Usually, though, we write $\chi$ as a combination of shape functions
# over the reference domain:
# $$
#   \chi(X) = \sum_{i=1}^m x_i N^e_i(X)
# $$
# where $x_i$ are given nodal locations.  In this case, the Jacobian of
# the map is
# $$
#   \chi'(X) = \sum_{i=1}^m x_i {N^e_i}'(X).
# $$
# When the same shape functions are used for defining the domain mapping
# and the interpolation of solution fields on the domain, we say we are
# using an *iso-parametric* mapping.

function isoparametric!(shapes, xnodal, x, J)
    shapes(x)
    mul!(x, xnodal, shapes.N,   1.0, 0.0)
    mul!(J, xnodal, shapes.dN', 1.0, 0.0)
    x
end

##
# It is useful to also define an isoparametric quadrature rule

struct IsoMappedRule{T <: QuadratureRule,
                     S <: ShapeFuns,
                     M <: AbstractMatrix} <: QuadratureRule
    base_rule :: T                # Reference domain rule
    shapes    :: S                # Shapes
    xnodal    :: M                # Nodal points
    map_grads :: Bool             # Flag whether to map gradients
    J         :: Matrix{Float64}  # Space for Jacobian/LU
    ipiv      :: Vector{Int}      # Pivots for Jacobian LU
end

IsoMappedRule(base_rule, shapes, xnodal, map_grads=false) =
    IsoMappedRule(base_rule, shapes, xnodal, map_grads,
                  zeros(quad_dim(base_rule), quad_dim(base_rule)),
                  zeros(Int, quad_dim(base_rule)))

r2s_gradients!(q :: IsoMappedRule, g) = r2s_gradients!(q.J, q.ipiv, g)

quad_npoints(q :: IsoMappedRule) = quad_npoints(q.base_rule)
quad_dim(q :: IsoMappedRule) = quad_dim(q.base_rule)

function quad_point(q :: IsoMappedRule, i)
    x = quad_point(q.base_rule, i)
    isoparametric!(q.shapes, q.xnodal, x, q.J)
    r2s_factor!(q.J, q.ipiv)
    if q.map_grads r2s_gradients!(q, q.shapes) end
    x
end

function quad_weight(q :: IsoMappedRule, i; map = false)
    if map quad_point(q, i) end
    quad_weight(q.base_rule, i) * r2s_det(q.J, q.ipiv)
end

