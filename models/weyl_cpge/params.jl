using SuperLattice.Util.Physics, SuperLattice.Util
using Statistics
αx, αy, αz, β, id = gammaMatrices(3; rep=:Dirac)

# model parameters
t=1.0
gamma=0.8
M=2
L=30
sizes=[L,L,L]
W=0.0
######## random function
rf = SuperLattice.Util.rand_gen(prod(sizes); rf=randn, offset=true) # random gaussian disorder that offset to 0
########


OBC=false
#tbc_theta = [0.5, 0.5, 0.5] # in unit of 2pi
# if not set, random twisted boundary condition set here.
tbc_theta = rand(3)
println("theta=$(tbc_theta)")



# KPM
NC=512
Ntilde = NC*6
NR=1

