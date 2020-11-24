using SuperLattice
using LinearAlgebra, Arpack, SparseArrays
using SuperLattice.Util.Physics


uc = UnitCell(;d=3)
addAtom(uc, :a, 2, [0, 0, 0]) # spin up / down

# potential term - disorder and skyrmion texture
Vr = function(r1, r2, v)
    # x and y are the same
    v .= W * rf() * σ0 
end
addPot(uc, :a, Vr)

# hopping term
addHopExt(uc, :a, :a, [1, 0, 0], 
          0.5im * t * σx; hc=true)
addHopExt(uc, :a, :a, [0, 1, 0], 
          0.5im * t * σy; hc=true)
addHopExt(uc, :a, :a, [0, 0, 1], 
          0.5im * t * σz; hc=true)

# set up lattice
ltc = ExplicitLattice(d=3,
                      pv=[[1., 0, 0],
                          [0., 1, 0],
                          [0., 0, 1]],
                      sizes=sizes,
                      OBC=OBC,
                      z_skip=false);
uc_symb = addUC(ltc, uc);
# add magnetic field
@time addMagneticField(ltc; B=B)

populateUC(ltc, uc_symb);

SuperLattice.refresh_none(ltc)
if isdefined(Main, :Hsp_gen)
    println("Hsp_gen defined, reusing")
    SuperLattice.pull_anonymous_funcs!(Hsp_gen, ltc)
else
    println("Hsp_gen not defined, creating new")
    Hsp_gen = get_operator_gen(ltc);
end
