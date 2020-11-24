"""
Define atom in the unit cell
atomSymb: a symbol (:a, :b, :c etc...)
dof: number of band/level/internal DOF: integer (default = 1)
idx_range: start/end of index when collecting atoms of the same unitcell together
xloc: ArrayTs: distance from the lattice point
Optional:
may use dof_cum_prev to determine idx_range rather than explicitly write it.
"""
struct Atom{Ts} 
    atomSymb::Symbol
    dof::Int64
    idx_range::UnitRange{Int64} # length = 2
    xloc::Array{Ts, 1}

    function Atom{Ts}(atomSymb::Symbol, dof::Integer, idx_range::UnitRange{T} where {T<:Integer}, xloc::Array) where {Ts}
        dof = Int64(dof)
        # check
        @assert dof > 0
        #        @assert (idx_range[2]-idx_range[1])==(dof-1) ##TODO repair this check
        new{Ts}(atomSymb::Symbol, dof, idx_range, convert(Array{Ts,1}, xloc))
    end
end
Atom{Ts}(atomSymb::Symbol, dof::Integer, idx_limits::Array{T,1} where {T<:Integer}, xloc::Array) where {Ts} =Atom{Ts}(atomSymb, dof, Int64(idx_limits[1]):Int64(idx_limits[2]), xloc)
Atom{Ts}(atomSymb::Symbol, dof::Integer, xloc::Array; dof_cum_prev::Integer=0) where {Ts} = Atom{Ts}(atomSymb, Int64(dof), [dof_cum_prev+1, dof_cum_prev+dof], xloc)
# if no name is given, assign based on xloc
function Atom{Ts}(dof::Integer, xloc::Array) where {Ts}
    println("Deprecated. Please assign specific name.")
    Atom{Ts}(Symbol(xloc), Int64(dof), xloc)
end


## Print functions
function Base.show(atom::Atom)
    println(typeof(atom), ", Name=", atom.atomSymb, ", DOF=", atom.dof, ", Loc=", atom.xloc)
end
Base.show(io::IO, a::Atom) = Base.show(a::Atom)
