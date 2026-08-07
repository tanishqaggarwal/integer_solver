"""Experiment 2: also close a7930 and a29539 via their free inputs."""
import sys, time
import dlib as L
import fwd

P = L.P
v = L.load('../best/new_instance_partial_39026.json')
for t in L.definer:
    v[t] = 0
fwd.forward(v)
v[9118] = 0; v[8731] = 0; v[1329] = 0; v[10903] = 0
fwd.forward(v)

for it in range(6):
    # a7930: 9367949*(x24548 - x25442) - x7927 ; x7927 = p*x11052
    v[24548] = v[25442]
    v[11052] = 0
    # a29539: 12846437*(x14853 - x1308) - x29967 ; x29967 = p*x30163
    v[14853] = v[1308]
    v[30163] = 0
    fwd.forward(v)
    av, nz, f = fwd.report(v, f'iter{it}', show=8)
    if len(nz) <= 2:
        break
L.save(v, 'D_state1.json')
