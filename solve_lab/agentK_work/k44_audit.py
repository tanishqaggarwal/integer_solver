#!/usr/bin/env python3
"""K44: re-run the three load-bearing closure-derived results under the GLOBAL forward guard.

Everything in this directory that reads a value out of a closure was produced without the
per-slot guard.  These three are the ones other conclusions rest on."""
import sys, os, json, random
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import fold as FD
from cascadep import CascadeP, NV, P
from k26_drive import FORBID
from k43_forward import pin, close, e2s, comp

C2 = CascadeP()
vc = json.load(open(K + '/varclass2.json'))
h, ls, ob, wr = vc['handles'], vc['leafsel'], vc['otherbools'], vc['wires']
dv = [u for u in range(NV) if u not in set(C2.E.free)]
S = FD.SHIFT


def run(on, obval, guard, pre=None):
    seed = {u: 0 for u in h}
    for u in ls: seed[u] = 1 if u in on else 0
    if obval is None:
        order = h + ls + dv + wr + ob
    else:
        for u in ob: seed[u] = obval
        order = h + ls + ob + dv + wr
    if pre:
        for u, val in pre.items(): seed[u] = val
        order = h + list(pre) + order[len(h):]
    v, _ = C2.close(seed, order, forbid=FORBID, pin=pin if guard else None)
    return v


print('=== A. does the dead-input gate stay off?  (k38 re-run, guarded) ===')
for obv, lab in ((0, 'other bools = 0'), (1, 'other bools = 1'), (None, 'derived')):
    v = run(set(), obv, True)
    bad, nz = C2.score(v)
    print('  all selectors OFF, %-16s root gate x15298 = %d   failing eqs = %d  atoms = %d'
          % (lab, v[15298], bad, len(nz)))

print('\n=== B. are the 900 non-leaf booleans inert?  (k35 re-run, guarded) ===')
for ON in ([0, 1, 3], [3, 10]):
    on = set(e2s[e] for e in ON)
    vals = []
    for obv in (0, 1, None):
        v = run(on, obv, True)
        vals.append(((v[12186] + S) % P, v[16742], (v[14853] + S) % P, v[24908]))
    print('  ON %-10s identical across bools=0/1/derived: %s' % (ON, vals[0] == vals[1] == vals[2]))

print('\n=== C. what does a wrong slot value cost, in EQUATIONS?  (k41 re-run, guarded) ===')
random.seed(5)
on = set(e2s[e] for e in (0, 1, 3))
v0 = run(on, 0, True); b0, nz0 = C2.score(v0)
print('  baseline ON={e0,e1,e3}: failing eqs = %d  atoms = %d' % (b0, len(nz0)))
for w, nm in ((12186, 'root A.x'), (16742, 'root A.y'), (14853, 'root B.x'), (24908, 'root B.y')):
    bad = random.randrange(2, P - 2)
    v = run(on, 0, True, pre={w: bad})
    b, nz = C2.score(v)
    print('  x%-6d %-9s kept=%-5s failing eqs = %-4d (delta %+d)  atoms = %d'
          % (w, nm, v[w] == bad, b, b - b0, len(nz)))
