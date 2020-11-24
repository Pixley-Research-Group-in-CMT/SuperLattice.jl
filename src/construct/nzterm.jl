## Data structure to record each non zero term

struct NZterm{Tv} 
    atomSymb_a::Symbol
    atomSymb_c::Symbol
    uc_c::Array{Int64, 1} # unit cell coordinate for creation operator
    t0::Symbol # symbol for look up in f_table

    function NZterm{Tv}(atomSymb_a::Union{Symbol, String}, atomSymb_c::Union{Symbol, String},uc_c::Array{T,1} where {T<:Integer}, t0::Symbol) where {Tv}
        uc_c = Int64.(uc_c)
        new{Tv}(Symbol(atomSymb_a), Symbol(atomSymb_c), uc_c, t0)
    end
end

## Display functions
function Base.show(a::NZterm)
    all_fields = fieldnames(typeof(a))
    for f in all_fields
        if !(occursin("Symb", String(f)))
            print(f, " = ", getfield(a, f), "; ")
        end
    end
end
Base.show(io::IO, a::NZterm) = Base.show(a::NZterm)
