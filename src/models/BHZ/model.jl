using KPMjulia, KPMjulia.Hamiltonians.Interface
using LinearAlgebra, Arpack, SparseArrays
using KPMjulia.Util.Physics


uc = UnitCell(;d=2)
addAtom(uc, :a, 2, [0,0]) # spin up / down



# potential term - disorder and skyrmion texture
Vr = function(r, r2, v)
    # x and y are the same
    v .= W * sum(cos.(Q .* r + φ)) * σ0 + (M-2) * σz
end
addPot(uc, :a, Vr)
# hopping term (only external)
addHopExt(uc, :a, :a, [-1, 0], 
          0.5 * t * (- 1im  * σx + σz); hc=true)
addHopExt(uc, :a, :a, [0, -1], 
          0.5 * t * (- 1im  * σy + σz); hc=true)

# set up lattice
pv = [
      [cos(ltc_y_rot_angle), sin(ltc_y_rot_angle)],
      [0.,1]
     ]

ltc = ExplicitLattice(d=2,pv=pv,sizes=sizes,OBC=OBC,z_skip=false);
uc_symb = addUC(ltc, uc);
# add magnetic field
# @time addMagneticField(ltc; B=[0,0,0.1])

populateUC(ltc, uc_symb);

Interface.refresh_none(ltc)
Hsp_gen = get_operator_gen(ltc);
