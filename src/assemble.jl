# --
# clear!(assembler) -- Clears things out
# assemble_add!(assembler, emat, ids) -- add matrix

using SparseArrays

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

clear!(A :: SparseMatrixCSC) = (A.nzval[:] .= 0.0)

function assemble_add!(A :: SparseMatrixCSC, emat, ids)
    nids = length(ids)
    for j = 1:nids
        if ids[j] >= 0
            
            # Start and end of col j in A and emat
            k, kn = A.colptr[ids[j]], A.colptr[ids[j]+1]-1
            i = 1
            
            # Merge pass (skip any not in graph)
            while k <= kn && i <= nids
                krow = A.rowval[k]
                irow = ids[i]
                if irow < krow # Includes case irow <= 0 (BC)
                    i += 1
                elseif krow < irow
                    k += 1
                else
                    A.nzval[k] += emat[i,j]
                    k += 1
                    i += 1
                end
            end
            
        end
    end
end
