using SuperLattice.Util, SuperLattice.Util.Physics, SuperLattice
using KPM
using LinearAlgebra
using DelimitedFiles, Random

Random.seed!(123)

include("params.jl") # this is the only part to change every time

include("model.jl")


println("create sparse matrix")
All_ops = to_sparse_matrix(Hsp_gen; Ops=[:H], da=[1.0,0.0,0],db=[0,1.0,0.0],tbc_theta=tbc_theta);
Hsp = All_ops[1];
println("max epsilon: $(maximum(abs.(Hsp-Hsp'))) (before explicit hermitianize)")
Hsp = (Hsp + Hsp') / 2 #
println("max epsilon: $(maximum(abs.(Hsp-Hsp')))")
println("size of Hsp: $(Hsp.m) x $(Hsp.n)")

All_ops = to_sparse_matrix(Hsp_gen; Ops=[:Ja], da=[1.0,0.0,0],db=[0,1.0,0.0],tbc_theta=tbc_theta);
Jx = All_ops[1];
All_ops = to_sparse_matrix(Hsp_gen; Ops=[:Jb], da=[1.0,0.0,0],db=[0,1.0,0.0],tbc_theta=tbc_theta);
Jy = All_ops[1];

a, H_norm = Util.normalizeH(Hsp);
println("a = $a")
@time mu = KPM.kpm_1d(H_norm, NC, 10, H_norm.n);
@time E, rho = KPM.dos(mu, a; N_tilde=Ntilde, E_range=[-3.0, 3.0]); 
writedlm("rho_cpge.dat", [E real(rho) imag(rho)])


@time μ2Dxx = KPM.kpm_2d(H_norm, Jx, Jx, NC,NR,H_norm.n; arr_size=206)
Evals = collect(-6.0:0.001:6.0)
@time dσxxE = KPM.d_dc_cond(μ2Dxx, a, Evals)
writedlm("DC-long_cpge.dat", [Evals real(dσxxE) imag(dσxxE)])
