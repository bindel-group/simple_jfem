include("../src/assemble.jl")

using Test

function build_tridiag(A, n)
    emat = [1.0 -1.0 ; -1.0 1.0]
    ids = [-1; 1:n-1; -1]

    clear!(A)
    for e = 1:n
        assemble_add!(A, emat, view(ids,e:e+1))
    end
    
    A
end

T5 = [ 2.0 -1.0  0.0  0.0  0.0 ;
      -1.0  2.0 -1.0  0.0  0.0 ;
       0.0 -1.0  2.0 -1.0  0.0 ;
       0.0  0.0 -1.0  2.0 -1.0 ;
       0.0  0.0  0.0 -1.0  2.0 ]

T5dense = build_tridiag(zeros(5,5), 6)
@test T5dense == T5

T5coo = build_tridiag(COOAssembler(24, 5, 5), 6)
T5csc = to_csc(T5coo)
@test T5csc == T5

T5csc = build_tridiag(T5csc, 6)
@test T5csc == T5
