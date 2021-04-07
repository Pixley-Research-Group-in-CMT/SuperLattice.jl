# Structure: Unit cell
#
# This section is not computational expensive, and
# should keep correctness and easiness of use at priority.
using LinearAlgebra, Logging
export UnitCell, addAtom, addPot, addHop, addHopInt, addHopExt

include("atom.jl")
include("nzterm.jl")

"""
The model in each unit cell. Only use the exposed 
functions to build model.
"""
struct UnitCell{Ts <: Real, Tv <: Complex}  
    # Types of value:
    # Ts is for spatial arguments. Always real. 
    # Tv is for values. Complex for most occasion.
    # Ts and Tv are expected to have matching 
    # precision

    # general parameters:
    a::Ts                                   # unit length
    d::Int64                                  # dimension of space
    atoms::Dict{Symbol, Atom{Ts}}   # info about each atom

    ## terms: 
    t_ex::Dict{NamedTuple{(:atom_a, :atom_c, :uc_c), 
                          Tuple{Symbol,Symbol,Array{Int64}}}, 
               NZterm{Tv}}    

    # table for functions
    f_table::Dict{Symbol, Function}

    # keep track of terms added. Init as 0
    counts::Array{Int64, 1} # first term NNZ (by terms), second term total dof


    # Initiate empty
    function UnitCell{Ts, Tv}(;a::Real=1.0, d::Integer=2) where {Ts, Tv}
        d = Int64(d)

        atoms = Dict{Symbol, Atom{Ts}}()
        t_ex = Dict{NamedTuple{(:atom_a, :atom_c, :uc_c), 
                               Tuple{Symbol,Symbol,Array{Int64}}}, 
                    NZterm{Tv}}()
        f_table = Dict{Symbol, Function}()
        new{Ts, Tv}(Ts(a), Int64(d), atoms, t_ex, f_table, [0, 0])
    end
end
UnitCell(args...; kwargs...) = UnitCell{Float64, ComplexF64}(args...; kwargs...)


get_xloc(uc::UnitCell, atom_symb::Symbol) = get_xloc(uc.atoms[atom_symb])
get_xloc!(uc::UnitCell, atom_symb::Symbol, xloc::Union{Array, SubArray}) = get_xloc!(uc.atoms[atom_symb], xloc)
add_xloc!(uc::UnitCell, atom_symb::Symbol, xloc::Union{Array, SubArray}) = add_xloc!(uc.atoms[atom_symb], xloc)


## Print functions
Base.show(io::IO, uc::UnitCell) = Base.show(uc)
function Base.show(uc::UnitCell)
    println("Unit Cell:")
    println("Spatial Dimension = ", uc.d)
    println("**Atoms: ")
    for x in keys(uc.atoms)
        println((x,), "==>", uc.atoms[x])
    end
    println("**Non Zero Terms (NZterm):")
    for x in keys(uc.t_ex)
        println(values(x), " ==> ", uc.t_ex[x])
    end
    print("---------------")
end


"""
add atom to uc
"""
function addAtom(uc::UnitCell{Ts, Tv}, atomSymb::Symbol, dof::Int64, xloc) where {Ts, Tv}
    @assert (uc.d == length(xloc)) "spatial dimension missmatch"    
    @assert (!(atomSymb in keys(uc.atoms))) "naming conflict"

    for v in values(uc.atoms)
        if v.xloc == xloc
            @warn "existing atom at this loc"
        end
    end

    atom = Atom{Ts}(atomSymb, dof, xloc; dof_cum_prev=get_dof(uc))
    uc.atoms[atom.atomSymb] = atom
    uc.counts[2] += dof # add to total dof
    println("Adding Atom: \n",atom)
end


""" 
add a potential term in uc, at atomSymb.
"""
function addPot(uc::UnitCell{Ts, Tv}, atomSymb::Symbol, v0::Function; prefix="p", hc=true) where {Ts, Tv}

    nargs_v0 = nargs_f(v0)
    if nargs_v0 == 0 # v0 has no arg
        t0 = function(r1, r2, v_out)
            v_out .= v0()
        end
    elseif nargs_v0 == 1 # v0 has one arg
        t0 = function(r1, r2, v_out)
            v_out .= v0(r1)
        end
    elseif nargs_v0 == 2 # v0 has two args
        t0 = function(r1, r2, v_out)
            v_out .= v0(r1, r2)
        end
    elseif nargs_v0 == 3 # v0 has 3 args
        t0 = v0
    else
        println("WARNING: addPot statement has bad v0. Skipping.")
        return
    end
    addTerm(uc, atomSymb, atomSymb, zeros(Int64, uc.d), t0; prefix=prefix, hc=hc)

end

function addPot(uc::UnitCell{Ts, Tv}, atomSymb::Symbol, v0::Array{T, 2} where {T<:Number}) where {Ts, Tv}
    v0_f = (xloc_a -> v0)
    return addPot(uc, atomSymb, v0_f)
end
function addPot(uc::UnitCell{Ts, Tv}, atomSymb::Symbol, v0::Number) where {Ts, Tv}
    atom_dof = get_dof(uc, atomSymb)
    if atom_dof == 1
        v0_mat = ones(Tv, 1,1) 
    else
        println("Warning: definition of potential V($(atomSymb)) = $(v0) is ambiguous. \n Taking V($(atomSymb)) = $(v0) I$(subscript(atom_dof)).")
        v0_mat = Array{Tv, 2}(I(atom_dof))
    end
    @.  v0_mat .*= v0
    return addPot(uc, atomSymb, v0_mat)
end


"""
Add term. Only "outbound", meaning `atom_a` (annihilation) is alwasy in `this`
unit cell, and `atom_c` (creation) is in unit cell `uc_c` away, in terms of
iloc (integer factor of primitive vectors). If `uc_c==0`, it is hopping to
another atom in the same unitcell. If `uc_c==0` and `atom_a==atom_c`, it is
a potential term.

`t0` is a function that defines the hopping based on the location of `atom_a`
and `atom_c`. `f` is not executed in the model building phase.

`hc` tells whether hermitian conjugate is automatically added.  By default
true. If selected to be false, hermitian conjugate term should be added
manually. When adding a potential term, hermitian conjugacy is enforced on `t0`


"""
function addTerm(
                uc::UnitCell{Ts, Tv}, 
                atomSymb_a::Symbol, 
                atomSymb_c::Symbol, 
                uc_c::Array{Int64,1}, 
                t0::Function;
                hc::Bool = true,
                prefix="he"
               ) where {Ts, Tv}
    @assert (atomSymb_a in keys(uc.atoms)) "atom_a not found - atomSymb_a=$(atomSymb_a)"
    @assert (atomSymb_c in keys(uc.atoms)) "atom_c not found - atomSymb_c=$(atomSymb_c)"

    # fix specific situations of pot, int and ext
    ispot = false
    if all(uc_c .== 0)
        if atomSymb_a == atomSymb_c
            println("term $(atomSymb_a) -> $(atomSymb_c) at $(uc_c) is potential")
            if prefix == "he"
                prefix = "p" 
            elseif prefix != "p"
                prefix = join(["p", prefix])
            end
            println("prefix reset as $(prefix)")
            ispot = true
            
            #if hc
            #    @assert check_hc(uc, t0, get_dof(uc, atomSymb_a)) "potential term $(atomSymb_a) is not Hermitian as expected"
            #end
        else
            println("term $(atomSymb_a) -> $(atomSymb_c) at $(uc_c) is hopping in same U.C.")
            if prefix == "he"
                prefix = "hi" 
            elseif prefix != "hi"
                prefix = join(["hi", prefix])
            end
            println("prefix reset as $(prefix)")
        end
    else
        if !startswith(prefix, "he")
            prefix = join(["he", prefix])
        end
    end

    @assert (length(uc_c)==uc.d) "Unitcell dimension $(d) mismatch with hopping distance uc_c dimension $(length(uc_c))"


    t0_symb = addFunction(uc, t0; prefix=prefix)
    uc.t_ex[(atom_a=atomSymb_a,atom_c=atomSymb_c,uc_c=uc_c)] = NZterm{Tv}(atomSymb_a, atomSymb_c, uc_c, t0_symb)
    uc.counts[1] += 1 # update counter for terms
    if (hc & !ispot)
        t0_adj_symb = addFunction(uc, t0; symb=Symbol(join([string(t0_symb), "hc"])), adjoint_out=true)
        uc.t_ex[(atom_a=atomSymb_c,atom_c=atomSymb_a,uc_c=-uc_c)] = NZterm{Tv}(atomSymb_c, atomSymb_a, -uc_c, t0_adj_symb)
        uc.counts[1] += 1 # update counter for terms
    end
end

function addTerm(uc::UnitCell{Ts, Tv},
                 atomSymb_a::Symbol,
                 atomSymb_c::Symbol,
                 uc_c::Array{Int64,1},
                 t0::Array{T,2} where {T<:Number};
                 kwargs...
                ) where {Ts, Tv}
    t0_f = function(xloc_a, xloc_c, t0_out)
        @. t0_out = t0
    end
    return addTerm(uc, atomSymb_a, atomSymb_c, uc_c, t0_f;kwargs...)
end

function addTerm(uc::UnitCell{Ts, Tv},
                 atomSymb_a::Symbol,
                 atomSymb_c::Symbol,
                 uc_c::Array{Int64,1},
                 t0::Number;
                 kwargs...
                ) where {Ts, Tv}
    atom_dof_a = get_dof(uc, atomSymb_a)
    atom_dof_c = get_dof(uc, atomSymb_c)
    if atom_dof_a == atom_dof_c == 1
        t0_mat = ones(Tv, 1,1) 
    else
        @assert (atom_dof_a == atom_dof_c) "Error: definition of potential t($(atomSymb_a)->$(atomSymb_c)) = $(t0) is ambiguous. Use a $(atomSymb_c) x $(atomSymb_a) matrix instead"
        t0_mat = Array{Tv, 2}(I(atom_dof_a))
    end
    @. t0_mat .*= t0
    return addTerm(uc, atomSymb_a, atomSymb_c, uc_c, t0_mat; kwargs...)
end

# aliases
addHopInt(uc, atomSymb_a, atomSymb_c, t0; kwargs...) = addTerm(uc, atomSymb_a, atomSymb_c, zeros(Int64, uc.d), t0; kwargs...)
addHopExt(args...; kwargs...) = addTerm(args...; kwargs...)
addHop(args...; kwargs...) = addTerm(args...; kwargs...)



# The nonzero entries starting (annihilating in) each uc 
# counted by the number of terms
get_term_count(uc::UnitCell) = uc.counts[1]


function get_nnz_per_uc(uc::UnitCell)
    nnz = 0

    nnz += sum(Int64[get_dof(uc, x.atom_a) * get_dof(uc, x.atom_c) for x in keys(uc.t_ex)])

    return nnz
end

function get_dof(uc::UnitCell, atom_symb::Symbol)
    return uc.atoms[atom_symb].dof
end

# the DOF of all atoms
get_dof(uc::UnitCell) = uc.counts[2]

get_idx_range(uc::UnitCell, atom_symb::Symbol) = uc.atoms[atom_symb].idx_range

"""
add a new symbol - function pair to function table
return the name of the function.
"""
function addFunction(uc::UnitCell, f::Function; symb::Symbol=:UNDEF, prefix="f", adjoint_out=false)
    if symb == :UNDEF
        for i in 1:(length(uc.f_table)+1)
            symb = Symbol(join([prefix, i]))
            if !haskey(uc.f_table, symb)
                return addFunction(uc, f; symb=symb, prefix=prefix)
            end
        end
    end
    
    output_modifier = identity
    if adjoint_out
        output_modifier = adjoint
    end
    # require each function to be single method
    @assert (length(methods(f))==1) "function $symb has more than one method."
    nargs = nargs_f(f)
    @assert (2<=nargs<=3) "function should either be f(x,y) -> v or f!(x, y, v)"
    if (nargs == 2)
        println("Improvement note: rewrite function $symb into inplace form function(xloc_a, xloc_c, v_out) to reduce memory allocation.")
        f! = function(xloc_a, xloc_c, v_out)
            v_out .= output_modifier(f(xloc_a, xloc_c))
        end
        uc.f_table[symb] = f!
    elseif (nargs == 3)
        f! = function(xloc_a, xloc_c, v_out)
            f(xloc_a, xloc_c, output_modifier(v_out))
        end
        uc.f_table[symb] = f!
    end

    return symb
end

function check_hc(uc::UnitCell, f::Function, dof::Int64; test_N::Int64=10)
    nargs = nargs_f(f)
    @assert (2<=nargs<=3) "function should either be f(x,y) -> v or f!(x, y, v)"
    x = zeros(Float64, uc.d)
    y = zeros(Float64, uc.d)
    v = zeros(Float64, dof, dof)
    if (nargs == 2)
        for i = 1:test_N
            x .= rand(Float64, uc.d)
            y .= rand(Float64, uc.d)
            if !(ishermitian(f(x, y)))
                return false
            end
        end
    else
        for i = 1:test_N
            x .= rand(Float64, uc.d)
            y .= rand(Float64, uc.d)
            f(x, y, v)
            if !ishermitian(v)
                return false
            end
        end
    end
    return true
end
