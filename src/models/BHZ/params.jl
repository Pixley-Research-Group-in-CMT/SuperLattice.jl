using KPMjulia.Util.Physics, KPMjulia.Util
using Statistics

# model parameters
t=1.0
L=89
sizes=[L,L]
W=0.00
M=3.00
Q=2*pi*34/89

OBC=false
tbc_theta = rand(rng, 2)
φ=rand(rng, 2) *2*pi

NC=1024 # 2^14
NR=13
Ntilde=NC * 3

ltc_y_rot_angle=0.0
a_fixed=4.0
