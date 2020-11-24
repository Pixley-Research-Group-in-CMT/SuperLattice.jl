using KPMjulia, KPMjulia.Hamiltonians.Interface
using LinearAlgebra, Arpack, SparseArrays
using KPMjulia.Util.Physics


uc = UnitCell(;d=2)
addAtom(uc, :a, 1, [0,0]) # spin up / down

function QPμ(r_a, r_c)
    J = 0.0im
    for (W, Q, φ) in zip(W_QP_HOP, Q_QP_HOP, phi_QP_HOP)
        J += W * sum(
                     cos.(Q * (r_a + r_c) / 2 + φ(r_a, r_c))
                    )
    end
    return J
end

# + pi/2 is for a convention that makes all parameters real.
PiFluxAx(r) = pi/2 + pi/2
PiFluxAy(r) = - (-1)^r[1] * pi/2 + pi/2


addHopExt(uc, :a, :a, [1, 0], 
          (ra, rc) -> (- (t + QPμ(ra, rc)) * exp.(1im * PiFluxAx(rc)));
          hc=true)
addHopExt(uc, :a, :a, [0, 1], 
          (ra, rc) -> (- (t + QPμ(ra, rc)) * exp.(1im * PiFluxAy(rc)));
          hc=true)

# set up lattice
ltc = ExplicitLattice(d=2,pv=[[1.,0],[0.,1]],sizes=sizes,OBC=OBC,z_skip=false);
uc_symb = addUC(ltc, uc);
# add magnetic field
# @time addMagneticField(ltc; B=[0,0,0.1])

populateUC(ltc, uc_symb);

Interface.refresh_none(ltc)
Hsp_gen = get_operator_gen(ltc);
