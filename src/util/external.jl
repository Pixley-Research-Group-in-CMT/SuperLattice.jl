# utility functions exposed to user under SuperLattice.Util
using Random
export rand_gen
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
