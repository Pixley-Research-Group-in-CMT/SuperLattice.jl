#################
# Abstract terms for unit cell definition
# #NOTE to myself: may want to move to another file
abstract type AbstractTerms{Tv <: Union{Real, Complex}} end # terms in hamiltonian

# Hopping = spatial different
#   internal = same UC
#   external = different UC
# Potential = spatial same (may cross bands)

abstract type AbstractHopping{Tv} <: AbstractTerms{Tv} end
abstract type AbstractHoppingInternal{Tv} <: AbstractHopping{Tv} end
abstract type AbstractHoppingExternal{Tv} <: AbstractHopping{Tv} end
abstract type AbstractPotential{Tv} <: AbstractTerms{Tv} end

ishop(t::AbstractHopping) = true
ishop(p::AbstractPotential) = false

ispot(p::AbstractHopping) = false
ispot(p::AbstractPotential) = true

issameUC(t::AbstractHoppingInternal) = true
issameUC(p::AbstractPotential) = true
issameUC(t::AbstractHoppingExternal) = false
################



abstract type AbstractAtom{Ts <: Real} end
get_xloc(atom::AbstractAtom) = atom.xloc
function get_xloc!(atom::AbstractAtom, xloc::Union{Array, SubArray})
    xloc .= atom.xloc
end
function add_xloc!(atom::AbstractAtom, xloc::Union{Array, SubArray})
    xloc .+= atom.xloc
end





#################
# Abstract definition for geometry
#

abstract type AbstractLattice{Ts <: Real} end



abstract type Operator{Ts <: Real, Tv <: Union{Real, Complex}} end
