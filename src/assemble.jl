using SparseArrays

##
# #
# Each element in a finite element discretization consists of
#
# - A domain $\Omega_e$ for the $e$th element, and
# - Local shape functions $N^e_1, \ldots, N^e_m$, which are often
#   Lagrange functions for interpolation at some set of nodes
#   in $\Omega_e$.
#
# Each local shape function on the domain $\Omega_e$ is the
# restriction of some global shape function on the whole
# domain $\Omega$.  That is, we have global shape functions
# $$
#   N_{j}(x) = \sum_{j = \iota(j',e)} N^e_{j'}(x),
# $$
# where $\iota(j,e)$ denotes the mapping from the local shape function
# index for element $e$ to the corresponding global shape function
# index.  We only ever compute explicitly with the local functions
# $N^e_j$; the global functions are implicit.
#
# *Assembly* is the process of reconstructing a quantity defined in
# terms of global shape functions from the contributions of the
# individual elements and their local shape functions.  For example,
# to compute
# $$
#   F_i = \int_{\Omega} f(x) N_i(x) \, dx,
# $$
# we rewrite the integral as
# $$
#   F_i = \sum_{i = \iota(i',e)} \int_{\Omega_e} f(x) N^e_{i'}(x) \, dx.
# $$
# In code, this is separated into two pieces:
#
# - Compute element contributions $\int_{\Omega_e} f(x) N^e_{i'}(x) \, dx$.
#   This is the responsibility of the element implementation.
# - Sum contributions into the global position $i$ corresponding to
#   the element-local index $i'$.  This is managed by an assembly loop.
#
# The concept of an "assembly loop" is central to finite element
# methods, but it is not unique to this setting.  For example, circuit
# simulators similarly construct system matrices (conductance,
# capacitance, etc) via the contributions of circuit elements
# (resistors, capacitors, inductors, and so forth).
#
# We have two types of assembly loops that we care about: those that
# involve pairs of shape functions and result in matrices, and those
# that explicitly involve only a single shape function and result in
# vectors.
#
# Our assemblers all implement two methods
#
# - `clear!(assembler)` -- Clears things out
# - `assemble_add!(assembler, econtrib, ids)` -- add the element
#   contribution `econtrib` at the locations indicated by `ids`.  Any
#   zero or negative indices are dropped.
#
# ## Filtered loops
#
# We will sometimes also want to discard some element contributions
# that correspond to interactions with shape functions associated with
# known boundary values (for example).  We also handle this filtering
# work as part of our assembly process.  Because we do this a lot,
# we define an iterator over valid `(i, id)` pairs, where `i` is the
# index in the element numbering system and `id` is the corresponding
# reduced index in the global system (with Dirichlet BC indices skipped).

struct IdIterator{T <: AbstractVector} ids :: T end

function Base.iterate(iter :: IdIterator, state=1)
    while state <= length(iter.ids) && iter.ids[state] <= 0
        state += 1
    end
    state > length(iter.ids) ? nothing : ((state, iter.ids[state]), state+1)
end

##
# ## Vector assembly

clear!(v :: Vector) = (v[:] .= 0)

function assemble_add!(v :: Vector, evec, ids)
    for (i, id) in IdIterator(ids)
        v[id] += evec[i]
    end
end

##
# ## Dense matrix assembly

clear!(A :: Matrix) = (A[:] .= 0)

function assemble_add!(A :: Matrix, emat, ids)
    for (j, idj) in IdIterator(ids)
        for (i, idi) in IdIterator(ids)
            A[idi,idj] += emat[i,j]
        end
    end
end

##
# ## Coordinate form assembly
#
# A coordinate form matrix (COO) is just a list of `(i,j,A[i,j])`
# tuples (which we will store in parallel arrays).  We will follow the
# convention that entries with duplicate row/column indices are summed
# in the final matrix.  We preallocate space for a certain number of
# entries, and keep a counter `nentries` for how much of that
# preallocation is used.  We can reallocate if needed, but of course
# it is better not to do so.

mutable struct COOAssembler
    I :: Vector{Integer}  # Row ids
    J :: Vector{Integer}  # Column ids
    V :: Vector{Float64}  # Values
    nentries :: Integer   # Number of entries saved
    m :: Integer          # Rows in matrix
    n :: Integer          # Columns in matrix
end

COOAssembler(nalloc :: Integer, m, n=0) =
    COOAssembler(zeros(Integer, nalloc),
                 zeros(Integer, nalloc),
                 zeros(nalloc), 0, m, n > 0 ? n : m)

##
# Clearing the coordinate form assembler doesn't require filling any
# space with zeros -- we just reset the `nentries` counter.

clear!(assembler :: COOAssembler) = (assembler.nentries = 0)

##
# The `ensure_capacity!` function ensures that we have capacity for
# `ncontribs` more entries.  If we do not have capacity in the pre-allocated
# space, we resize the arrays to either double the preallocated capacity
# or to accommodate the additional contributions, whichever is more.

function ensure_capacity!(assembler :: COOAssembler, ncontribs)
    nold = length(assembler.V)
    if assembler.nentries + ncontribs > nold
        nnew = max(assembler.nentries + ncontribs, 2*nold)
        resize!(assembler.I, nnew)
        resize!(assembler.J, nnew)
        resize!(assembler.V, nnew)
    end
end

##
# The `add_entry!` function adds a single entry, assuming that capacity
# has already been ensured.

function add_entry!(assembler :: COOAssembler, idi, idj, entry)
    assembler.nentries += 1
    assembler.I[assembler.nentries] = idi
    assembler.J[assembler.nentries] = idj
    assembler.V[assembler.nentries] = entry
end

##
# Finally, the `assembler_add!` function ensures that we have enough
# capacity, and then adds all the entries from the element matrix.

function assemble_add!(assembler :: COOAssembler, emat, ids)
    ensure_capacity!(assembler, length(ids)^2)
    for (j, idj) in IdIterator(ids)
        for (i, idi) in IdIterator(ids)
            add_entry!(assembler, idi, idj, emat[i,j])
        end
    end
end

##
# The Julia `SparseArrays` package has a built-in function already to
# convert a coordinate form matrix in parallel arrays into a compressed
# sparse column representation.

to_csc(a :: COOAssembler) =
    sparse(view(a.I, 1:a.nentries),
           view(a.J, 1:a.nentries),
           view(a.V, 1:a.nentries), a.m, a.n)

##
# ## Compressed sparse column reassembly
#
# The `CSCAssembler` keeps the data for assembling a compressed sparse
# column matrix.  In addition to the storage for the matrix that we are
# assembling, we keep some auxiliary scratch storage for aggregating
# contributions to one column.  This scratch storage is initialized to
# zeros, and should stay all zeros outside the `CSCAssembler` routines.

struct CSCAssembler{Tv,Ti}
    A :: SparseMatrixCSC{Tv,Ti}  # CSC storage structure
    Aj :: Vector{Tv}             # Scratch vector
end

CSCAssembler(A :: SparseMatrixCSC{Tv,Ti}) where {Tv,Ti} =
    CSCAssembler(A, zeros(Tv, A.m))

function clear!(a :: CSCAssembler)
    a.A.nzval[:] .= 0.0
    a.Aj[:] .= 0.0
end

function assemble_add!(a :: CSCAssembler, emat, ids)
    nids = length(ids)
    A, Aj = a.A, a.Aj
    for (j, idj) in IdIterator(ids)

        # Populate dense column scratch
        for (i, idi) in IdIterator(ids)
            Aj[idi] += emat[i,j]
        end

        # Extract from dense column
        k1, kn = A.colptr[ids[j]], A.colptr[ids[j]+1]-1
        for k = k1:kn
            A.nzval[k] += Aj[A.rowval[k]]
        end

        # Clear dense columns scratch
        for (i, idi) in IdIterator(ids)
            Aj[idi] = 0.0
        end
    end
end
