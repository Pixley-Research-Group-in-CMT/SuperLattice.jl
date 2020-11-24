using KPMjulia, KPMjulia.Hamiltonians.Interface
using LinearAlgebra, Arpack, SparseArrays
using KPMjulia.Util.Physics

αx, αy, αz, β, id = gammaMatrices(3; rep=:Dirac)

uc = UnitCell(;d=3)
addAtom(uc, :a, 4, [0,0,0]) # 4 dof each point

# potential term - disorder
Vr = function(x, y, v)
    v .= W * rf() * id + (m0 + 3 * m2) * β
end
addPot(uc, :a, Vr)
# hopping term (only external)
addHopExt(uc, :a, :a, [1, 0, 0], 
          0.5im * t * αx - 0.5 * m2 * β; hc=true)
addHopExt(uc, :a, :a, [0, 1, 0], 
          0.5im * t * αy - 0.5 * m2 * β; hc=true)
addHopExt(uc, :a, :a, [0, 0, 1], 
          0.5im * t * αz - 0.5 * m2 * β; hc=true)

# set up lattice
ltc = ExplicitLattice(d=3,pv=[[1.,0,0],[0.,1,0],[0.,0,1]],sizes=sizes,OBC=OBC);
uc_symb = addUC(ltc, uc);
# add magnetic field
# @time addMagneticField(ltc; B=[0,0,0.1])
populateUC(ltc, uc_symb);
Interface.refresh_none(ltc)
if isdefined(Main, :Hsp_gen)
    println("Hsp_gen defined, reusing")
    Interface.pull_anonymous_funcs!(Hsp_gen, ltc)
else
    println("Hsp_gen not defined, creating new")
    Hsp_gen = get_operator_gen(ltc);
end
println("Finishing ltc.")
