using SuperLattice.Util.Physics, SuperLattice.Util
using Statistics

# KPM
NC=512
Ntilde = NC*6
NR=15

# model parameters
t=1.0
L=34
sizes=[L,L]
W_QP_HOP=[-0.3, 1.5]
Q_QP_HOP=[pi, 13/34*2*pi]
phi_QP_HOP = [[0.5*pi, 0.5*pi], rand(2)*pi]

OBC=true
tbc_theta = [0.0, 0.0]

num2f(x) = replace(string(round(x, digits=3)), "-"=>"m")

W_str = join(num2f.(W_QP_HOP), "-")
Q_str = join(num2f.(Q_QP_HOP), "-")
fn_params = "L$(L)_W$(W_str)_Q$(Q_str)"
