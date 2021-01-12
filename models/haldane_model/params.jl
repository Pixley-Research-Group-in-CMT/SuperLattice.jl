using SuperLattice.Util.Physics, SuperLattice.Util
using Statistics

# KPM
NC=4096
Ntilde = NC*6
NR=15

# model parameters
t1=1.0
t2=0.5
M=0.2
L=233
sizes=[L,L]

OBC=false
tbc_theta=rand(rng, 2)

fn_params = "L$(L)_M$(M)_t1$(t1)_t2$(t2)"
