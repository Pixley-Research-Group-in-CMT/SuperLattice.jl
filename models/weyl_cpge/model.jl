using SuperLattice
using LinearAlgebra, Arpack, SparseArrays
using SuperLattice.Util.Physics


println("define model")
uc = UnitCell(;d=3)
addAtom(uc, :a, 2, [0,0,0]) # 2 band model

# potential term - disorder
Vr = function(x, y, v)
    v .= W * rf() * σ0 + M * σz
end
addPot(uc, :a, Vr)
# hopping term (only external)
addHop(uc, :a, :a, [1, 0, 0], 
          -0.5im * t * σx - 0.5 * t * σz ; hc=true)
addHop(uc, :a, :a, [0, 1, 0], 
          -0.5im * t * σy - 0.5 * t * σz; hc=true)
addHop(uc, :a, :a, [0, 0, 1], 
          0.5im * gamma * σ0 - 0.5 * t * σz; hc=true)

# set up lattice
ltc = Lattice(d=3,pv=[[1.,0,0],[0.,1,0],[0.,0,1]],sizes=sizes,OBC=OBC);
uc_symb = addUC(ltc, uc);
# add magnetic field
# @time addMagneticField(ltc; B=[0,0,0.1])

populateUC(ltc, uc_symb);
SuperLattice.refresh_none(ltc)
Hsp_gen = get_operator_gen(ltc);
