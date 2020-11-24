mod_onebase(x, y) = mod(x-1, y)+1

function wrap_idx(iloc::Array{Int64, 1},
                  sizes::Array{Int64, 1};
                  z_skip=false)
    iloc_layer = sum(iloc)
    iloc_mod = mod_onebase.(iloc, sizes)
    iloc_mod[end] += z_skip * (iloc_layer - iloc_mod[1] - iloc_mod[2] - iloc_mod[end])
    return iloc_mod
end


function get_all_ind(sizes::Array{T,1} where {T<:Int64}, ext_dims::Int64)
    f(x) = collect(1:x)
    return collect(Iterators.product(map(f, sizes[1:ext_dims])...))[:]
end
get_all_ind(sizes::Array{T,1} where {T<:Int64}) = get_all_ind(sizes, length(sizes))


function ijk2I(ijk::Array{T,1}where{T<:Int64}, sizes::Array{T,1}where{T<:Int64}, ext_dims::Int64;
              z_skip::Bool=false)
    sizes = sizes[1:ext_dims]
    ijk = wrap_idx(ijk, sizes; z_skip=z_skip)
    #ijk = Int64.(mod_onebase.(ijk, sizes))
    mult = 1
    I = 0
    for (i, ind) in enumerate(ijk)
        I += (ind-1)*mult
        mult *= sizes[i]
    end
    return I+1
end
function I2ijk(I::Int64, sizes::Array{T,1} where {T<:Int64}, ext_dims::Int64; z_skip::Bool=false)
    I -= 1
    ijk = ones(Int64, ext_dims)
    for (i,l) in enumerate(sizes[1:ext_dims])
        ijk[i] = mod(I, l)+1
        I = div(I, l)
    end
    return ijk
end

I2ijk(I::Int64, sizes::Array{T,1}where{T<:Int64}) = I2ijk(I, sizes, length(sizes))

ijk2I(ijk::Array{T,1}where{T<:Int64}, sizes::Array{T,1} where {T<:Int64}) = ijk2I(ijk, sizes, length(sizes))
