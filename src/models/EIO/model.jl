using KPMjulia, KPMjulia.Hamiltonians.Interface
using DelimitedFiles



## LOAD DATA
include("atoms.jl")

fn_in = "EIO_symWAN_U$(U)_Rspace"

bands = 24
steps = bands ^ 2


pv = [[5.1350, 5.1350, 0.0],
      [0.0, 5.1350, 5.1350],
      [5.1350, 0.0, 5.1350]]

A = readdlm(fn_in);
tot_lines = length(@view A[:,1])

n_all = A[:,1:3]
unique_n = unique(n_all; dims=1)

band_n_all = A[:,4]
band_m_all = A[:,5]
unique_band = unique(band_m_all)

t = A[:,6] + 1im * A[:,7]


## BUILD UNITCELL
Z = 31+15
uc = UnitCell(; d=3)
ltc = ExplicitLattice(d=3,pv=pv,sizes=[L,L,Z],OBC=false,z_skip=true);
for (i, a_symb) in enumerate(atom_symbs)
    xloc = sum(atom_ilocs[i].*pv)
    addAtom(uc, a_symb, atom_dof[i], xloc)
end
visited = Set(Array[])
count = [0, 0]
for i in 1:Int64(tot_lines / steps)
    uc_to = Int64.(unique_n[i, :])
    if (-uc_to) in visited
        continue
    end
    push!(visited, uc_to)

    idx = ((i - 1) * steps + 1):(i*steps)
    band_n = band_n_all[idx]
    band_m = band_m_all[idx] # m is annihilation
    val = t[idx]
  
    for (atom_a, atom_c) in Iterators.product(atom_symbs, atom_symbs)
        #truncation
        xloc_a = Interface.get_xloc(uc, atom_a)
        xloc_c = Interface.get_xloc(uc, atom_c) + Interface.get_xloc(ltc, uc_to)
        if sqrt(sum((xloc_c-xloc_a).^2)) > truncate_trsd
            count[1] += 1
            continue
        else
            count[2] += 1
        end


        band_a = unique_band[orbit_atoms .== atom_a]
        band_c = unique_band[orbit_atoms .== atom_c]
        
        idx_single_term = map(nm -> ((nm[1] in band_c) & (nm[2] in band_a)), zip(band_n, band_m) )

        val_single_term = val[idx_single_term]
        val_mat_single_term = reshape(val_single_term, length(band_a), length(band_c))
        println("$atom_a -> $atom_c , $(uc_to), $(size(val_mat_single_term))")
        addHopExt(uc, atom_a, atom_c, uc_to, val_mat_single_term; hc=true)
    end

end


## BUILD LATTICE
uc_symb = addUC(ltc, uc);
#l0 = Int64(ceil((3*L-layer)/2))
#l1 = l0 + layer - 1
l0 = 2*L+2
l1 = l0 + layer - 1
@assert ((l1-1-1) < Z) "The boundary condition may go weird"

for ii in 1:L
    for jj in 1:L
        for kk in 1:Z
            ll = ii + jj + kk
            if (ll >= l0) & (ll <= l1)
                println("INCLUDE: $ii, $jj, $kk; layer=$ll;")
                populateUC(ltc, uc_symb; range=:SINGLE, loc=[ii,jj,kk])
            end
#            if (abs(ii-ll/3)+abs(jj-ll/3)+abs(kk-ll/3))> (0.5*L)
#                println("REMOVE: $ii, $jj, $kk; layer=$ll;")
#                continue
#            end
        end
    end
end

# add magnetic field
# addMagneticField(ltc; B=[0,0,0.1])
#@time populateUC(ltc, uc_symb);
Interface.refresh_none(ltc);


if isfile("Hgen.jld2")
    Hsp_gen = Interface.loadSMG("Hgen.jld2", ltc)
else
    @time Hsp_gen = get_operator_gen(ltc);
    Interface.nnz_compress!(Hsp_gen)
    Interface.nt_compress!(Hsp_gen)
    Interface.idx_compress!(Hsp_gen)
    Interface.saveSMG(Hsp_gen, "Hgen.jld2")
    # temp workaround: explicitly pull functions back. 
    Interface.pull_anonymous_funcs!(Hsp_gen, ltc)
end
