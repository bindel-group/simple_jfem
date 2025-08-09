#ldoc on
##
# # Quadrature rules
#
# Quadrature rules approximate integrals with formulas of the form
# $$
#   \int_{\Omega} f(x) \, d\Omega(x) \approx
#   \sum_{j=1}^p f(\xi_{j}) w_j
# $$
# where $\xi_j \in \Omega$ and $w_j \in \mathbb{R}$ are known as the
# quadrature nodes (or points) and weights, respectively.
#
# A good source of quadrature rules for various domains can be found
# in Stroud's book on _Approximate calculation of multiple integrals_
# (Prentice Hall, 1971).
# 
# ## Gaussian-Legendre quadrature rules
# 
# Gauss-Legendre quadrature rules (sometimes just called Gauss
# quadrature rules when the context is clear) are $p$-point rules on
# $[-1, 1]$ that are characterized by the fact that they are exact
# when $f$ is a polynomial of degree at most $2p-1$.
# 
# Gauss-Legendre nodes are zeros of Legendre polynomials, while the
# weights can be computed via an eigenvalue decomposition (using the
# Golub-Welsch algorithm).  However, we do not need very high-order
# quadrature rules, and so only provide nodes and weights for rules up
# to $p = 10$ (probably more than we need), which are tabulated in
# many places.  Because this is just a table lookup, we don't bother
# to include the code in the automated documentation.

function gauss_point(i, npts)
    gauss_pts = (
        # ... the points
#ldoc on
        # One point
        0.0,

        # Two points
        -0.577350269189626,
        0.577350269189626,

        # Three points
        -0.774596669241483,
        0.0,
        0.774596669241483,

        # Four points
        -0.861136311594053,
        -0.339981043584856,
        0.339981043584856,
        0.861136311594053,

        # Five points
        -0.906179845938664,
        -0.538469310105683,
        0.0,
        0.538469310105683,
        0.906179845938664,

        # Six points
        -0.932469514203152,
        -0.661209386466265,
        -0.238619186083197,
        0.238619186083197,
        0.661209386466265,
        0.932469514203152,

        # Seven points
        -0.949107912342759,
        -0.741531185599394,
        -0.405845151377397,
        0.0,
        0.405845151377397,
        0.741531185599394,
        0.949107912342759,

        # Eight points
        -0.960289856497536,
        -0.796666477413627,
        -0.525532409916329,
        -0.183434642495650,
        0.183434642495650,
        0.525532409916329,
        0.796666477413627,
        0.960289856497536,

        # Nine points
        -0.968160239507626,
        -0.836031107326636,
        -0.613371432700590,
        -0.324253423403809,
        0.0,
        0.324253423403809,
        0.613371432700590,
        0.836031107326636,
        0.968160239507626,

        # Ten points
        -0.973906528517172,
        -0.865063366688985,
        -0.679409568299024,
        -0.433395394129247,
        -0.148874338981631,
        0.148874338981631,
        0.433395394129247,
        0.679409568299024,
        0.865063366688985,
        0.973906528517172
#ldoc off
      )

    gauss_pts[( npts*(npts-1) )/2 + i]
end

function gauss_weight(i, npts)
    gauss_wts = (
        # ... and the weights
#ldoc off
        # One point
        2.0,

        # Two points
        1.0,
        1.0,

        # Three points
        0.555555555555556,
        0.888888888888889,
        0.555555555555556,

        # Four points
        0.347854845137454,
        0.652145154862546,
        0.652145154862546,
        0.347854845137454,

        # Five points
        0.236926885056189,
        0.478628670499366,
        0.568888888888889,
        0.478628670499366,
        0.236926885056189,

        # Six points
        0.171324492379170,
        0.360761573048139,
        0.467913934572691,
        0.467913934572691,
        0.360761573048139,
        0.171324492379170,

        # Seven points
        0.129484966168870,
        0.279705391489277,
        0.381830050505119,
        0.417959183673469,
        0.381830050505119,
        0.279705391489277,
        0.129484966168870,

        # Eight points
        0.101228536290376,
        0.222381034453374,
        0.313706645877887,
        0.362683783378362,
        0.362683783378362,
        0.313706645877887,
        0.222381034453374,
        0.101228536290376,

        # Nine points
        0.081274388361574,
        0.180648160694857,
        0.260610696402935,
        0.312347077040003,
        0.330239355001260,
        0.312347077040003,
        0.260610696402935,
        0.180648160694857,
        0.081274388361574,

        # Ten points
        0.066671344308688,
        0.149451349150581,
        0.219086362515982,
        0.269266719309996,
        0.295524224714753,
        0.295524224714753,
        0.269266719309996,
        0.219086362515982,
        0.149451349150581,
        0.066671344308688
#ldoc on        
    )

    gauss_wts[( npts*(npts-1) )/2 + i]
end

##
# ## Quadrature iterator interface
#
# The `QuadratureRule` abstract type provides a base type for all quadrature
# rules.  We assume that it provides methods
#
# - `quad_npoints(rule)`: Returns the number of points
# - `quad_point(rule, i)`: Returns the quadrature point $\xi_i$
# - `quad_weight(rule, i)`: Returns the quadrature weight $w_i$
#
# One can also provide a `quad_pointwt` that returns the point and weight
# as a pair (by default, this just calls the individual `quad_point` and
# `quad_weight` methods.
#
# For any quadrature rule, we overload the `Base.iterate` interface so that
# we can use it in `for` loops, writing expressions like
# ```{.julia}
# I = sum( f(xi) * wt for (xi, wt) in rule )
# ```
# to approximate the integral via the rule.

abstract type QuadratureRule end

quad_pointwt(q :: QuadratureRule, i) =
    (quad_point(q, i), quad_weight(q, i))

Base.iterate(q :: QuadratureRule, state=1) =
    state > quad_npoints(q) ? nothing : (quad_pointwt(q, state), state+1)

##
# ## 1D Gauss quadrature interfaces
# 
# Our `GaussRule1d` rule just provides an alternate interface to the
# `gauss_point` and `gauss_weight` functions defined earlier

struct GaussRule1d <: QuadratureRule
    npts :: Integer
end

quad_npoints(q :: GaussRule1d) = q.npts
quad_point(q :: GaussRule1d, i) = gauss_point(i, q.npts)
quad_weight(q :: GaussRule1d, i) = gauss_weight(i, q.npts)

##
# The `GaussRule1dv` rule returns the quadrature points in a vector of
# length 1 (rather than returning them as scalars)

struct GaussRule1dv <: QuadratureRule
    xi :: Vector{Float64}
    npts :: Integer
end

GaussRule1dv(npts :: Integer) = GaussRule1dv(zeros(1), npts)

quad_npoints(q :: GaussRule1dv) = q.npts

function quad_point(q :: GaussRule1dv, i)
    q.xi[1] = gauss_point(i, q.npts)
    q.xi
end

quad_weight(q :: GaussRule1dv, i) = gauss_weight(i, q.npts)

##
# ## Product Gauss rules
# 
# A 2D tensor product Gauss rule for the domain $[-1,1]^2$ involves a
# grid of `npts1`-by-`npts1` quadrature points with coordinates given
# by 1D Gauss quadrature rules.

struct GaussRule2d <: QuadratureRule
    xi :: Vector{Float64}
    npts1 :: Integer
 end

GaussRule2d(npts1 :: Integer) = GaussRule2d(zeros(2), npts1)

quad_npoints(q :: GaussRule2d) = q.npts1 * q.npts1

function quad_point(q :: GaussRule2d, i)
    ix, iy = ((i-1)%q.npts1)+1, div(i-1,q.npts1)+1
    q.xi[:] .= (gauss_point(ix, d), gauss_point(iy, d))
    q.xi
end

function quad_weight(q :: GaussRule2d, i)
    ix, iy = ((i-1)%q.npts1)+1, div(i-1,q.npts1)+1
    gauss_weight(ix, d) * gauss_weight(iy, d)
end

function quad_pointwt(q :: GaussRule2d, i)
    ix, iy = ((i-1)%q.npts1)+1, div(i-1,q.npts1)+1
    q.xi[:] .= (gauss_point(ix, q.npts1), gauss_point(iy, q.npts1))
    wt = gauss_weight(ix, q.npts1) * gauss_weight(iy, q.npts1)
    (q.xi, wt)
end

##
# ## Triangle mid-side rule
# For a triangle, a rule based on the three mid-side values is exact
# for every polynomial with total degree less than or equal to 2
# (which is enough for our purposes).  This is sometimes called the
# Hughes formula.

struct HughesRule2d <: QuadratureRule
    xi :: Vector{Float64}
end

HughesRule2d() = HughesRule2d(zeros(2))

quad_npoints(q :: HughesRule2d) = 3
quad_weight(q :: HughesRule2d, i) = 1.0/6

function quad_point(q :: HughesRule2d, i)
    pts = (0.5, 0.0,
           0.5, 0.5,
           0.0, 0.5)
    q.xi[:] .= pts[2*i-1:2*i]
end
