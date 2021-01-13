# utility functions exposed to user under SuperLattice.Util
using Random
export rand_gen, normalizeH

"""
returns a function that output an element of an array of N pre-generated random numbers.
Allow offsetting the average of the N numbers to 0
"""
function rand_gen(N::Integer; rank::Integer=0, rf::Function=randn, offset::Bool=false, rng=nothing)
    if !offset 
        println("WARNING: untested.  Please use native random functions for now. ")
    end
    if rng == nothing
        rng = Random.GLOBAL_RNG
    end

    rn_all = rf(rng, N)
    for p in 1:rank # when paralellized, replicate rank
        rn_all .= rf(rng, N)
    end

    if offset
        rn_mean = sum(rn_all) / length(rn_all)
        rn_all .-= rn_mean
    end

    ptr = [0]
    function rf_out()
        ptr[1] = mod(ptr[1], N) + 1
        return rn_all[ptr[1]]
    end
    return rf_out
end


"""
Normalize H. If requested, allow renormalizing it to fixed value.
"""
function normalizeH(H::SparseMatrixCSC{ComplexF64,Int64}; eps::Float64=0.1, setA::Float64=0.0)
    
    println("hermitian check: ")
    @assert abs(sum(H-H')) < 1e-16*sqrt(H.n)
    println("pass.")

	if setA==0
    es, vs = eigs(H;tol=0.001,maxiter=300)
#	println(es)
    Emax = maximum(abs.(es))
    Emin = -Emax
    a = (Emax - Emin)/(2 - eps)
	else
    a = setA
	end

    H = H/a
    return a, H

end
