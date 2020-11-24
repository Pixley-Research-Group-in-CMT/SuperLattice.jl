# struct lattice encodes where each uc is located
#
INF = 999999 # for periodic boundary condition. need to think about how to make it more sensible
export Lattice, addUC, populateUC, addMagneticField
struct Lattice{Ts} 
    d::Int64                              # dimension of space
    pv::Array{Array{Ts,1}}                # primitive vectors
    sizes::Array{Int64, 1}
    OBC::Array{Bool, 1}
    uc_table::Dict{Symbol, UnitCell}    # look up table for unit cell
    uc_map::Array{Symbol}               # map of UC symbols, where idx is loc of the uc
    counts::Array{Int64, 1}
    z_skip::Bool

    # Initiate empty
    function Lattice{Ts}(;d::Integer=2,sizes::Array{T,1} where {T<:Integer}=ones(Int64,2), OBC::Union{Array{Bool,1},Bool}=false,pv::Array{Array{Ts,1}}=[[0.,1.],[1.,0.]],z_skip=false) where Ts
        d = Int64(d)
        sizes = Int64.(sizes)
        if length(OBC)==1
            OBC = ones(Bool, d) * OBC
        end
        @assert (length(sizes) == d) "dimension mismatch: sizes parameter should have length $(d)"
        @assert (length(OBC) == d) "dimension mismatch: OBC parameter should have length $(d)"
        @assert (length(pv)== d ==length(pv[1])) "dimension mismatch: there should be d=$(d) primitive vectors with dimension $(d)"
        uc_table = Dict{Symbol, UnitCell}()
        uc_map = fill(:none, sizes...)
        new{Ts}(d, pv, sizes,OBC, uc_table, uc_map, [0,0], z_skip)
    end
end
Lattice(args...;kwargs...) = Lattice{Float64}(args...;kwargs...)

function addMagneticField(ltc::Lattice; B::Array{T, 1} where {T<:Real}=[0,0,0.0])
    for uc in values(ltc.uc_table)
        addMagneticField(uc; B=B)
    end
    return nothing
end
function addMagneticField(uc::UnitCell; B::Array{T, 1} where {T<:Real}=[0,0,0.0])
    peierls = gen_peierls(; B=B, d=uc.d)
    for (f_symb, f) in uc.f_table
        if startswith("p", string(f_symb))
            continue
        end
        function f_new(r1, r2, v)
            f(r1, r2, v)
            v .*= peierls(r1, r2)
        end
        uc.f_table[f_symb] = f_new
    end
    return nothing
end

function addUC(ltc::Lattice, uc::UnitCell)
    @assert (ltc.d==uc.d) "Lattice dimension $(ltc.d) mismatch with UC dimension $(uc.d)"
    # not given a name
    for i in 1:(length(ltc.uc_table)+1)
        uc_symb = Symbol(join(["UC", i]))
        if !haskey(ltc.uc_table, uc_symb)
            addUC(ltc, uc_symb, uc)
            return uc_symb
        end
        i += 1
    end
    println("Failed. Please specify a symbol to name the uc")
end


function addUC(ltc::Lattice, uc_symb::Symbol, uc::UnitCell)
    if (uc_symb => uc) in ltc.uc_table
        println("Already exist")
    else
        @assert !haskey(ltc.uc_table, uc_symb) "$(uc_symb) already representing an different existing UC"
        @assert !(uc in values(ltc.uc_table)) "The uc is already represented by $(filter(x->ltc.uc_table[x]==uc, keys(ltc.uc_table)))"
        ltc.uc_table[uc_symb] = uc
    end
    return uc_symb
end
# TODO: add renameUC

function populateUC(ltc::Lattice, uc_symb::Symbol, uc::UnitCell; range::Symbol=:ALL)
    addUC(ltc, uc_symb, uc)
    populateUC(ltc, uc_symb, range)
end

function populateUC(ltc, uc_symb::Symbol; range::Symbol=:ALL, loc::Array{T,1} where {T<:Integer}=[1,1])
    loc = Int64.(loc)
    @assert (uc_symb in keys(ltc.uc_table)) "UnitCell $(uc_symb) not found"
        
    if range==:ALL
        all_idx = get_all_ind(ltc)
        for idx in all_idx
            ltc.uc_map[idx...] = uc_symb
        end
        ltc.counts .+= length(all_idx) * ltc.uc_table[uc_symb].counts
    elseif range==:SINGLE
        ltc.uc_map[loc...] = uc_symb
        ltc.counts .+= ltc.uc_table[uc_symb].counts
    else
        println("range $(range) - WIP.")
    end
end


# get OBC info from lattice data structure
get_uc_symb(ltc::Lattice, 
            iloc::Array{T,1} where {T<:Integer}
            ) = get_uc_symb(ltc, iloc, ltc.OBC)

# explicit OBC info by Boolean
get_uc_symb(ltc::Lattice, 
            iloc::Array{T,1} where {T<:Integer},
            OBC::Union{Array{Bool,1},Bool}
            ) = get_uc_symb(ltc, 
                            iloc, 
                            ltc.sizes.+INF*(OBC),ltc.z_skip)
# explicit OBC info by number
function get_uc_symb(ltc::Lattice, 
                     iloc::Array{T,1} where {T<:Integer}, 
                     modulo::Union{Array{Int64, 1}, Int64} ,
                     z_skip::Bool=false)
    ### NOTE: z_skip is a temporary fix. Try to integrate into the code later. 
#    iloc_mod = mod.(iloc .- 1, modulo) .+ 1 # -1 and +1 for range 1 to modulo
    #iloc_layer = sum(iloc) # when z_skip is on, take final z from layer
    #iloc_mod = mod_onebase.(iloc, modulo)
    # if x, y changed, z up by z_skip
    # WARNING:: Turning on z_skip in 2D will screw up everything.
#    iloc_mod[end] += z_skip * (iloc_mod[1]!=iloc[1] | iloc_mod[2]!=iloc_mod[2])
#    iloc_mod[end] = mod_onebase(iloc_mod[end], modulo[end])
    #iloc_mod[end] += z_skip * (iloc_layer - iloc_mod[1] - iloc_mod[2] - iloc_mod[end]) # determine z from layer
    iloc_mod = wrap_idx(iloc, modulo; z_skip=z_skip)        


#    println(iloc, iloc_mod)
#@time    iloc_mod = Int64.(iloc_mod)
    
#    uc_symb = ltc.uc_map[iloc_mod...]
#@time    bool = all(iloc_mod .<= ltc.sizes) & all(iloc_mod .>= 1)
    bool = all(i -> 1 <= iloc_mod[i] <= ltc.sizes[i], eachindex(iloc_mod))
    iloc_mod .*= bool
    iloc_mod .+= !bool # if true, these two operations do not affect; if false, take it to 1,1,1
    ret = (bool, ltc.uc_map[iloc_mod...])

    #println("$(ret), $iloc_mod")
    return ret

#    if any(iloc_mod .> ltc.sizes) || any(iloc_mod .< 1)
##        println("location $(iloc_mod) out of bound")  ##INFO
#        return (false, :none)
#    end
#    
#    if isassigned(ltc.uc_map, iloc_mod...)
#        return (true, ltc.uc_map[iloc_mod...])
#    else
##        println("location $(iloc_mod) is Empty.") ##INFO
#        return (false, :none)
#    end
end

function get_uc(ltc::Lattice, uc_symb::Symbol)
    return ltc.uc_table[uc_symb]
end

get_all_uc_symb(ltc::Lattice) = keys(ltc.uc_table)
#get_uc_count(ltc::Lattice) = prod(ltc.sizes)
get_uc_count(ltc::Lattice) = sum(ltc.uc_map .!= :none)
get_all_ind(ltc::Lattice) = get_all_ind(ltc.sizes, ltc.d)

function get_xloc(ltc::Lattice, iloc::Array{Int64,1})
#    @assert (length(iloc)==ltc.d) "iloc should match dimension $(ltc.d)"
#    return sum(ltc.pv .* iloc)
    xloc = similar(ltc.pv[1])
    get_xloc!(ltc, iloc, xloc)
    return xloc
end
function get_xloc!(ltc::Lattice, iloc::Array{Int64,1}, xloc::Union{Array, SubArray})
    xloc .= 0
    add_xloc!(ltc, iloc, xloc)
end
function add_xloc!(ltc::Lattice, iloc::Array{Int64,1}, xloc::Union{Array, SubArray})
    for i in eachindex(ltc.pv)
        xloc .+= (iloc[i] * ltc.pv[i])
    end
end

get_term_count(ltc::Lattice) = ltc.counts[1]
get_dof(ltc::Lattice) = ltc.counts[2]


## Print functions
Base.show(io::IO, ltc::Lattice) = Base.show(ltc)
function Base.show(ltc::Lattice)
    println("Lattice:")
    println("Spatial Dimension = ", ltc.d)
    println("**UnitCells (uc_table):\n $(keys(ltc.uc_table))")
    println("**UnitCell Map (uc_map):")
    display(ltc.uc_map)
end


# consolidate f_table
function collect_f_table!(ltc, new_f_table)
    for (uc_symb, uc) in ltc.uc_table
        for (f_symb, f) in uc.f_table
            new_f_table[Symbol(join([string(uc_symb), string(f_symb)]))] = f
        end
    end
end

"""
Should not be explicitly invoked.

This function creates in ltc a unitcell for `:none`.
"""
function refresh_none(ltc)
    new_none = UnitCell(;d=ltc.d)
    for (uc_symb, uc) in ltc.uc_table
        if uc_symb == :none
            continue
        end
        for (atom_symb, atom) in uc.atoms
            if atom_symb in keys(new_none.atoms)
                println("$atom_symb repeat")
                continue
            end
            addAtom(new_none, atom_symb, atom.dof, atom.xloc)
        end
    end
    ltc.uc_table[:none] = new_none
end


# estimate number of nonzero terms 
function est_nnz(ltc::Lattice)
    # use maximum nnz per uc
    all_uc_symb = get_all_uc_symb(ltc)
    mx_nnz_per_uc = 0
    for uc_symb in all_uc_symb
        uc = get_uc(ltc, uc_symb)
        mx_nnz_per_uc = max(mx_nnz_per_uc, get_nnz_per_uc(uc))
    end

    # multiply by number of uc
    return mx_nnz_per_uc * get_uc_count(ltc)
end


# spy
import UnicodePlots.spy
function spy(ltc::Lattice)
    uc_map_num = similar(ltc.uc_map, Int64)
    for (i, uc_symb) in enumerate(keys(ltc.uc_table))
        println("$(uc_symb) -> $i")
        uc_map_num[ltc.uc_map .== uc_symb] .= i
    end
    println("none -> 0")
    uc_map_num[ltc.uc_map .== :none] .= 0
    if ltc.d == 2
        spy(uc_map_num)
    elseif ltc.d == 3
        for i in 1:ltc.sizes[3]
            println("layer $i:\n ----------------\n ")
            println(spy(uc_map_num[:,:,i], maxwidth=20,maxheight=20))
        end
    end
end


