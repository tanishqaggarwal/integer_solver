"""S11 step 35: targeted construction in frame 3, branch (1,1).

With x_7075 = 0:  a35759 = -x_29854 and a35761 = x_31864, so BOTH pairs close by
setting x_29854 = x_1329 = 0 and x_31864 = x_10903 = 0 -- a two-variable move the
greedy cannot see.  Then close a22229/a22230/a35762 through the detached
parameters, and the branch-activated gadgets a19088/a22233/a22235 through their
own handles.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame3 import DETACH, definer, ORDER, FREE, CHECKS, fwd, score, SSET
P = ad.P
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
v = list(base); v[2081] = 1; v[4287] = 1
fwd(v, rounds=10)
print(f'branch(1,1) frame3 raw {score(v)}')

# the four paired closures
v[29854] = 0; v[1329] = 0                     # a35758 = a35759 = 0
v[31864] = 0; v[10903] = 0                    # a35760 = a35761 = 0
v[642] = v[17325] * P                         # a35762 = 0
v[28730] = v[9413] * P                        # a22230 = 0
fwd(v, rounds=10)
v[4432] = v[19964] + v[28730]                 # a22231 = 0
fwd(v, rounds=10)
v[7068] = v[2099] + 7376877 * v[642]          # a22229 = 0
fwd(v, rounds=10)
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
print(f'after the paired closures: score {score(v)}  nonzero {nz}')

# the branch-activated gadgets, through their own handles / definers
for chk, gate, handles in [(22233, 22232, (6947, 26874)),
                           (22235, 22234, (24490, 33168)),
                           (19088, 19087, (950, 30095))]:
    if av[chk] == 0: continue
    ov = L.atom_out.get(gate)
    if ov is None: continue
    t = ov[1]
    tgt = T.solve_lin(chk, t, v)
    if tgt is None:
        print(f'  a{chk}: cannot solve for x_{t}'); continue
    vv = list(v); vv[t] = tgt
    got = None
    for h in handles:
        if h not in FREE: continue
        nv = T.solve_lin(gate, h, vv)
        if nv is not None:
            tr = list(v); tr[h] = nv
            fwd(tr, rounds=10)
            a2 = L.all_atom_values(tr)
            if a2[chk] == 0:
                v = tr; got = h; break
    av = L.all_atom_values(v)
    print(f'  a{chk}: {"closed via x_"+str(got) if got else "NOT closable"};  '
          f'score {score(v)}')
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
s = score(v)
print(f'\nFINAL score {s}  failing {L.NEQ-s}  nonzero {nz}')
if s > 39026:
    T.save(v, os.path.join(HERE, f'B11B_{s}.json'))
    print(f'  *** BEATS THE DELIVERABLE -- saved B11B_{s}.json')
