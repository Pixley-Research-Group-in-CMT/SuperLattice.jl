using KPMjulia.Util.Physics, KPMjulia.Util
using Statistics

# KPM
NC=4096
Ntilde = NC*6
NR=15

# model parameters
t=1.0
L=233
sizes=[L,L]
W_QP_POT=[0.3]
Q_QP_POT=[2pi * 89 / 233]
phi_QP_POT = rand(rng, 2) * pi

OBC=false
tbc_theta = rand(rng, 2) * pi

W_str = join(round.(W_QP_POT, digits=3), "-")
Q_str = join(round.(Q_QP_POT, digits=3), "-")
fn_params = "L$(L)_W$(W_str)_Q$(Q_str)"
