using KPMjulia.Util.Physics, KPMjulia.Util
using Statistics

# KPM
NC=512
Ntilde = NC*6
NR=1

# model parameters
t_nn=1.0
t_nnn_1=0.5 * t_nn
t_nnn_2=-0.2 * t_nn

# lattice parameters
L=100
sizes=[L,L]
OBC=true
## set tbc
#tbc_theta = [0.5, 0.5, 0.5] # in unit of 2pi
## if not set, random twisted boundary condition set here.
#tbc_theta = rand(2)
#for p in 1:(1*2*(prod(sizes)*NR)+1)
#    tbc_theta .= rand(2)
#end
## zero tbc
tbc_theta = rand(2) * 0

# modulation and disorder
W=0.0
######## random function
rf(x) = (rand(x) - 0.5) * 2#rand_gen(prod(sizes); rank=rank, rf=randn, offset=true) # random gaussian disorder that offset to 0
########
