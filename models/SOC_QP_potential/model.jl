using SuperLattice
using LinearAlgebra, Arpack, SparseArrays
using SuperLattice.Util.Physics


uc = UnitCell(;d=2)
addAtom(uc, :a, 2, [0,0]) # spin up / down

function QP(r)
    V = zeros(ComplexF64, 2, 2)
    for (W, Q, φ) in zip(W_QP_POT, Q_QP_POT, phi_QP_POT)
        V += W * σ0 * sum(
                          cos.(Q * r + φ)
                         )
    end
    return V
end
addPot(uc, :a, QP)
addHopExt(uc, :a, :a, [1, 0], 
          0.5im * t * σx;
          hc=true)
addHopExt(uc, :a, :a, [0, 1], 
          0.5im * t * σy;
          hc=true)

# set up lattice
ltc = Lattice(d=2,pv=[[1.,0],[0.,1]],sizes=sizes,OBC=OBC,z_skip=false);
uc_symb = addUC(ltc, uc);

# add magnetic field
# @time addMagneticField(ltc; B=[0,0,0.1])

populateUC(ltc, uc_symb);

SuperLattice.refresh_none(ltc)
Hsp_gen = get_operator_gen(ltc);
