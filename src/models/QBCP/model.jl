#https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.103.046811


using KPMjulia, KPMjulia.Hamiltonians.Interface
using LinearAlgebra, Arpack, SparseArrays
using KPMjulia.Util.Physics


uc = UnitCell(;d=2)
# use the blue (:a) as origin.
# set lattice size to 1
addAtom(uc, :a, 1, [0, 0]) 
addAtom(uc, :b, 1, [0.5, 0.5]) 


# potential term - disorder
Vr = function(x, y, v)
    # x and y are the same
    v .= W * rf()
end
addPot(uc, :a, Vr)

# hopping term 
addHopExt(uc, :a, :b, [0, 0], 
          - t_nn; hc=true)
addHopExt(uc, :a, :b, [0, -1], 
          - t_nn; hc=true)
addHopExt(uc, :a, :b, [-1, -1], 
          - t_nn; hc=true)
addHopExt(uc, :a, :b, [-1, 0], 
          - t_nn; hc=true)

addHopExt(uc, :a, :a, [0, 1], 
          - t_nnn_1; hc=true)
addHopExt(uc, :b, :b, [1, 0], 
          - t_nnn_1; hc=true)

addHopExt(uc, :b, :b, [0, 1], 
          - t_nnn_2; hc=true)
addHopExt(uc, :a, :a, [1, 0], 
          - t_nnn_2; hc=true)

# set up lattice
ltc = ExplicitLattice(d=2,pv=[[1.,0],[0.,1]],sizes=sizes,OBC=OBC,z_skip=false);
uc_symb = addUC(ltc, uc);
# add magnetic field
# @time addMagneticField(ltc; B=[0,0,0.1])

populateUC(ltc, uc_symb);

Interface.refresh_none(ltc)
Hsp_gen = get_operator_gen(ltc);
