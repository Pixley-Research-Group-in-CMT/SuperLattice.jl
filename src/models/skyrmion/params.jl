using KPMjulia.Util.Physics, KPMjulia.Util
using Statistics

# KPM
NC=512
Ntilde = NC*6
NR=1

# model parameters
t=1.0
L=100
sizes=[L,L]
W=0.3
######## random function
rf = rand_gen(prod(sizes); rank=rank, rf=randn, offset=true) # random gaussian disorder that offset to 0
########


OBC=true
#tbc_theta = [0.5, 0.5, 0.5] # in unit of 2pi
# if not set, random twisted boundary condition set here.
tbc_theta = rand(2)
for p in 1:(1*2*(prod(sizes)*NR)+1)
    tbc_theta .= rand(2)
end
tbc_theta = rand(2) * 0

# skyrmion profile
r_sk=10
Delta=0.67
m=1 # Vorticity. 1 for single skyrmion
gamma = pi/2 # 0 - Hedgehog; pi/2 - Vortex
γ=gamma


