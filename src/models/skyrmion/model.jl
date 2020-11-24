using SuperLattice
using LinearAlgebra, Arpack, SparseArrays
using SuperLattice.Util.Physics


uc = UnitCell(;d=2)
addAtom(uc, :a, 2, [0,0]) # spin up / down

θ(r) = 2 * pi - 4 * atan(exp(norm(r)/r_sk))
φxy(r) = angle(r[1] + 1im * r[2]) # azimuthal angle

mag(r) = [
          sin(θ(r)) * cos(m * φxy(r) + γ),
          sin(θ(r)) * sin(m * φxy(r) + γ),
          cos(θ(r))
         ]


# potential term - disorder and skyrmion texture
Vr = function(x, y, v)
    # x and y are the same
    v .= W * rf() * σz  - Delta * σ3dot(mag(x - [L/2, L/2]))
end
addPot(uc, :a, Vr)
# hopping term (only external)
addHopExt(uc, :a, :a, [1, 0], 
          - t * σ0; hc=true)
addHopExt(uc, :a, :a, [0, 1], 
          - t * σ0; hc=true)

# set up lattice
ltc = ExplicitLattice(d=2,pv=[[1.,0],[0.,1]],sizes=sizes,OBC=OBC,z_skip=false);
uc_symb = addUC(ltc, uc);
# add magnetic field
# @time addMagneticField(ltc; B=[0,0,0.1])

populateUC(ltc, uc_symb);

SuperLattice.refresh_none(ltc)
Hsp_gen = get_operator_gen(ltc);
