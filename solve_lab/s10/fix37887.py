"""S11 step 28: can a37887 be repaired?

In frame 3 x_28730 is free and costs exactly one equation: a37887.  Freeing A1
takes the twelve from 7 failing to 6, so if a37887 can be re-satisfied the total
is 6  ->  39,027.  a37887 is a single quadratic check over many variables; look
for any variable that can absorb it.
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
print(f'frame 3 base {score(base)}')
v = list(base)
v[28730] = v[9413] * P                       # A1 = 0
fwd(v, rounds=8)
av = L.all_atom_values(v)
print(f'A1 = 0: a22230 {av[22230]}, a22231 {"0" if av[22231]==0 else "NONZERO"}')
nz = [a for a in range(L.NA) if av[a]]
print(f'  nonzero {nz}  failing {len(L.failing_eqs(av))}  score {score(v)}')
print(f'  a37887 = {"ZERO" if av[37887]==0 else "nonzero"}')

vs = sorted(set(L.avars[37887]))
print(f'\na37887 has {len(vs)} variables; {sum(1 for w in vs if w in FREE)} are free')
ok = []
for w in vs:
    tgt = T.solve_lin(37887, w, v)
    if tgt is None or tgt == v[w]: continue
    cands = []
    if w in FREE: cands.append((w, tgt, 'direct'))
    else:
        d = definer.get(w)
        if d is not None:
            vv = list(v); vv[w] = tgt
            for u in sorted(set(L.avars[d])):
                if u == w or u not in FREE: continue
                nv = T.solve_lin(d, u, vv)
                if nv is not None: cands.append((u, nv, f'via x_{w}'))
    for u, nv, how in cands:
        tr = list(v); tr[u] = nv
        fwd(tr, rounds=8)
        at = L.all_atom_values(tr)
        nzt = [a for a in range(L.NA) if at[a]]
        s = score(tr)
        ok.append((s, u, how, nzt))
        if s > 39026:
            T.save(tr, os.path.join(HERE, f'WIN_{s}.json'))
            print(f'  *** x_{u} ({how}) -> SCORE {s}  nonzero {nzt}   SAVED')
ok.sort(reverse=True)
print(f'\n{len(ok)} repair attempts; best:')
for s, u, how, nzt in ok[:10]:
    print(f'  score {s:>6}  x_{u:<7} {how:<12} nonzero {nzt}')
