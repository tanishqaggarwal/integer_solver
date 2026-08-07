"""Experiment 1: forward-eval, then repair the four cheap modular checks directly."""
import sys, time
import dlib as L
import fwd

P = L.P
v = L.load('../best/new_instance_partial_39026.json')
for t in L.definer:
    v[t] = 0
fwd.forward(v)
av, nz, f = fwd.report(v, 'base', show=0)

# Step 1: zero x_9118, x_8731 (free) and their handles x_1329, x_10903
for u in (9118, 8731, 1329, 10903):
    print('x_%d free=%s val_len=%d' % (u, u in L.freeset, len(str(v[u]))))
v[9118] = 0
v[8731] = 0
v[1329] = 0
v[10903] = 0
fwd.forward(v)
fwd.report(v, 'after zero 9118/8731', show=20)
