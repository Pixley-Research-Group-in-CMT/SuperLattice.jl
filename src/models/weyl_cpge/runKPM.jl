using KPMjulia.Util, KPMjulia.Util.Physics, KPMjulia.Hamiltonians.Interface
using KPMjulia.KPM
using LinearAlgebra
using DelimitedFiles, Random

# parallel and distributed
n = parse(Int64,ARGS[1])
r = 0
Random.seed!(n)

include("params.jl") # this is the only part to change every time

println("theta=$(tbc_theta)")

println("define model")
include("model.jl")


println("create sparse matrix")
All_ops = Interface.to_sparse_matrix(Hsp_gen; Ops=[:H], da=[1.0,0.0,0],db=[0,1.0,0.0],tbc_theta=tbc_theta);
Hsp = All_ops[1];
println("max epsilon: $(maximum(abs.(Hsp-Hsp'))) (before explicit hermitianize)")
Hsp = (Hsp + Hsp') / 2 #
println("max epsilon: $(maximum(abs.(Hsp-Hsp')))")
println("size of Hsp: $(Hsp.m) x $(Hsp.n)")

All_ops = Interface.to_sparse_matrix(Hsp_gen; Ops=[:Ja], da=[1.0,0.0,0],db=[0,1.0,0.0],tbc_theta=tbc_theta);
Jx = All_ops[1];
All_ops = Interface.to_sparse_matrix(Hsp_gen; Ops=[:Jb], da=[1.0,0.0,0],db=[0,1.0,0.0],tbc_theta=tbc_theta);
Jy = All_ops[1];

a, H_norm = KPMjulia.Util.normalizeH(Hsp);
println("a = $a")
@time mu = KPMjulia.KPM.DOS_KPM(H_norm, NC, 10, H_norm.n);
@time E, rho = KPMjulia.KPM.computeDOSslow(mu, a, Ntilde;Erange=[-3.0, 3.0]); 
writedlm("rho_n$(n)_r$(r).dat", [E real(rho) imag(rho)])

#    @time μ2Dxy = KPMjulia.KPM.KPM_2D_FAST(H_norm, Jx, Jy, NC,NR,H_norm.n; arr_size=206)
#    @time Evals, dσxyE = KPMjulia.KPM.computeOCIntegrand(μ2Dxy, a, Ntilde; Erange=[-3.0,3.0],NC=0)
#    writedlm("DC-hall_n$(n)_r$(r).dat", [Evals real(dσxyE) imag(dσxyE)])


@time μ2Dxx = KPMjulia.KPM.KPM_2D_FAST(H_norm, Jx, Jx, NC,NR,H_norm.n; arr_size=206)
@time Evals, dσxxE = KPMjulia.KPM.computeOCIntegrand(μ2Dxx, a, Ntilde; Erange=[-3.0,3.0],NC=0)
writedlm("DC-long_n$(n)_r$(r).dat", [Evals real(dσxxE) imag(dσxxE)])
