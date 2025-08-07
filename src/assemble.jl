# --
# clear!(assembler) -- Clears things out
# assemble_add!(assembler, emat, ids) -- add matrix

using SparseArrays

# -- Vector assembly

clear!(v :: Vector) = (v[:] .= 0)

function assemble_add!(v :: Vector, evec, ids)
    nids = length(ids)
    for i = 1:nids
        if ids[i] > 0
            v[ids[i]] += evec[i]
        end
    end
end

# -- Dense matrix assembly

clear!(A :: Matrix) = (A[:] .= 0)

function assemble_add!(A :: Matrix, emat, ids)
    nids = length(ids)
    for j = 1:nids
        for i = 1:nids
            if ids[i] > 0 && ids[j] > 0
                A[ids[i],ids[j]] += emat[i,j]
            end
        end
    end
end

# -- Coordinate form assembly

mutable struct COOAssembler
    I :: Vector{Integer}
    J :: Vector{Integer}
    V :: Vector{Float64}
    nentries :: Integer
    m :: Integer
    n :: Integer
end

COOAssembler(nalloc :: Integer, m, n=0) =
    COOAssembler(zeros(Integer, nalloc),
                 zeros(Integer, nalloc),
                 zeros(nalloc), 0, m, n > 0 ? n : m)

clear!(assembler :: COOAssembler) = (assembler.nentries = 0)

function assemble_add!(assembler :: COOAssembler, emat, ids)

    # Allocate more space if needed
    nold = length(assembler.V)
    if assembler.nentries + length(ids)^2 > nold
        nnew = max(n + length(ids)^2, 2*nold)
        resize!(assembler.I, nnew)
        resize!(assembler.J, nnew)
        resize!(assembler.V, nnew)
    end

    # Add to the list
    for j = 1:length(ids)
        for i = 1:length(ids)
            if ids[i] > 0 && ids[j] > 0
                assembler.nentries += 1
                assembler.I[assembler.nentries] = ids[i]
                assembler.J[assembler.nentries] = ids[j]
                assembler.V[assembler.nentries] = emat[i,j]
            end
        end
    end

end

to_csc(a :: COOAssembler) = sparse(view(a.I, 1:a.nentries),
                                   view(a.J, 1:a.nentries),
                                   view(a.V, 1:a.nentries), a.m, a.n)

# -- Assemble into exisiting CSC matrix

struct CSCAssembler{Tv,Ti}
    A :: SparseMatrixCSC{Tv,Ti}
    Aj :: Vector{Tv}
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
    for j = 1:nids
        if ids[j] >= 0

            # Populate dense column scratch
            for i = 1:nids
                if ids[i] >= 0
                    Aj[ids[i]] += emat[i,j]
                end
            end

            # Extract from dense column
            k1, kn = A.colptr[ids[j]], A.colptr[ids[j]+1]-1
            rows = view(A.rowval, k1:kn)
            vals = view(A.nzval,  k1:kn)
            vals[:] .+= view(Aj,rows)

            # Clear dense columns scratch
            for i = 1:nids
                if ids[i] >= 0
                    Aj[ids[i]] = 0.0
                end
            end

        end
    end
end
