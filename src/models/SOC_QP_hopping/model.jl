using KPMjulia, KPMjulia.Hamiltonians.Interface
using LinearAlgebra, Arpack, SparseArrays
using KPMjulia.Util.Physics


uc = UnitCell(;d=2)
addAtom(uc, :a, 2, [0,0]) # spin up / down

function QPμ(r_a, r_c)
    J = 0.0im
    for (W, Q) in zip(W_QP_HOP, Q_QP_HOP)
        J += W * sum(
                     cos.(Q * (r_a + r_c) / 2 + φ)
                    )
    end
    return J
end

addHopExt(uc, :a, :a, [1, 0], 
          (ra, rc) -> (1im * (t + QPμ(ra, rc)) * σx); 
          hc=true)
addHopExt(uc, :a, :a, [0, 1], 
          (ra, rc) -> (1im * (t + QPμ(ra, rc)) * σy);
          hc=true)

# set up lattice
ltc = ExplicitLattice(d=2,pv=[[1.,0],[0.,1]],sizes=sizes,OBC=OBC,z_skip=false);
uc_symb = addUC(ltc, uc);
# add magnetic field
# @time addMagneticField(ltc; B=[0,0,0.1])

populateUC(ltc, uc_symb);

Interface.refresh_none(ltc)
Hsp_gen = get_operator_gen(ltc);
