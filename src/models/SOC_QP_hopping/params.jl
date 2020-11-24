using KPMjulia.Util.Physics, KPMjulia.Util
using Statistics

# KPM
NC=512
Ntilde = NC*6
NR=15

# model parameters
t=1.0
L=100
sizes=[L,L]
W_QP_HOP=[0.3]
Q_QP_HOP=[pi]

OBC=true
tbc_theta = [0.0, 0.0]
φ = [0.0, 0.0]

W_str = join(round.(W_QP_HOP, digits=3), "-")
Q_str = join(round.(Q_QP_HOP, digits=3), "-")
fn_params = "L$(L)_W$(W_str)_Q$(Q_str)"
