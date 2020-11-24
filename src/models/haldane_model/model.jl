####
# Tutorial ref: https://topocondmat.org/w4_haldane/haldane_model.html
# Original ref: https://doi.org/10.1103/PhysRevLett.61.2015
# Real space convention: in the unit of distance between nearest A-B sites
# 
#      A(0,1) -- B(0,1)
#     /            \
#   B(-1,1)        A(1,0)
#     \            /
#      A(0,0) -- B(0,0)
#     /            \
#   B(-1,0)        A(1,-1)

using KPMjulia, KPMjulia.Hamiltonians.Interface
using LinearAlgebra, Arpack, SparseArrays
using KPMjulia.Util.Physics


uc = UnitCell(;d=2)
addAtom(uc, :a, 1, [0,0]) # site A
addAtom(uc, :b, 1, [1,0]) # site B

# vanilla graphene model
addHop(uc, :a, :b, [0, 0], t1; hc=true)
addHop(uc, :a, :b, [-1, 0], t1; hc=true)
addHop(uc, :a, :b, [-1, 1], t1; hc=true)

# break symmetry
addPot(uc, :a, M)
addPot(uc, :b, -M)

# NNN term
addHop(uc, :a, :a, [1, 0], 1im * t2; hc=true)
addHop(uc, :a, :a, [0, 1], -1im * t2; hc=true)
addHop(uc, :a, :a, [1, -1], -1im * t2; hc=true)

addHop(uc, :b, :b, [-1, 0], 1im * t2; hc=true)
addHop(uc, :b, :b, [0, 1], 1im * t2; hc=true)
addHop(uc, :b, :b, [-1, 1], -1im * t2; hc=true)


# set up lattice
primitive_vectors = [
                     [1.5, sqrt(3)/2],
                     [0, sqrt(3)]
                    ]
ltc = ExplicitLattice(d=2,
                      pv=primitive_vectors,
                      sizes=sizes,
                      OBC=OBC,
                      z_skip=false); # z_skip is a temporary API. Set to false unless you know what it is. 
uc_symb = addUC(ltc, uc);

# add magnetic field
# @time addMagneticField(ltc; B=[0,0,0.1])

populateUC(ltc, uc_symb);

Interface.refresh_none(ltc)
Hsp_gen = get_operator_gen(ltc);
