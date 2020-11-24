using SuperLattice.Util.Physics, SuperLattice.Util
using Statistics
αx, αy, αz, β, id = gammaMatrices(3; rep=:Dirac)

# model parameters
t=1.0
m2=0.5
mratio=-1.87
m0 = m2 * mratio
L=50
sizes=[L,L,L]
W=0.6
######## random function
rf = rand_gen(prod(sizes); rf=randn, offset=true, rng=rng) # random gaussian disorder that offset to 0
########


OBC=false
tbc_theta = [0.5, 0.5, 0.5] # in unit of 2pi



# KPM
NC=1024
Ntilde = NC*2
NR=7
arr_size=224

a_fixed=5.0
