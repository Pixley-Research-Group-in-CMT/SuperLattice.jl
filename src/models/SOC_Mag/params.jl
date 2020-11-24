using KPMjulia.Util.Physics, KPMjulia.Util
using Statistics

# KPM
NC=1024
Ntilde = NC*6
NR=15
a = 4.0

# model parameters
t=1.0
L=100
sizes=[L, L, L]
W=0.8
Bz=0.1
B=[0, 0, Bz]
######## random function
rf = rand_gen(prod(sizes); rf=randn, offset=true, rng=rng) # random gaussian disorder that offset to 0
# Need two occurance per site
########


OBC=false
#tbc_theta = [0.5, 0.5, 0.5] # in unit of 2pi
# if not set, random twisted boundary condition set here.
tbc_theta = rand(3)

a_fixed=5.0
