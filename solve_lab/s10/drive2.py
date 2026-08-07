"""S11 step 33: cheapest live driver of x_19964.

The pair move needs x_19964 to shift by -d.  Scan every free input in its
structural ancestor cone: does it move x_19964, and what does it break outside
the eight?
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame3 import DETACH, definer, ORDER, FREE, CHECKS, fwd, score, SSET
P = ad.P
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
cone, st = set(), [19964]
while st:
    t = st.pop()
    if t in cone: continue
    cone.add(t)
    a = definer.get(t)
    if a is None: continue
    for w in L.avars[a]:
        if w != t: st.append(w)
fr = sorted(w for w in cone if w in FREE)
print(f'ancestor cone of x_19964: {len(cone)} vars, {len(fr)} free inputs', flush=True)
t0 = time.time()
hits = []
for u in fr:
    v = list(base); v[u] = v[u] + 1
    fwd(v, rounds=8)
    dv = v[19964] - base[19964]
    if dv == 0: continue
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a] and a not in SSET]
    eqs = set()
    for a in nz: eqs |= set(L.atom2eq[a])
    hits.append((len(eqs), u, dv, nz))
    print(f'  x_{u:<7} moves x_19964 by {dv};  outside-eight {nz} '
          f'({len(eqs)} eqs)', flush=True)
hits.sort()
print(f'\n{len(hits)} live drivers ({time.time()-t0:.0f}s)')
if hits:
    print(f'cheapest: x_{hits[0][1]} at {hits[0][0]} equations')
    # try the full pair move with the cheapest driver
    d = base[9413] * P - base[28730]
    for cost, u, dv, nz in hits[:4]:
        if d % dv: 
            print(f'  x_{u}: step not integral ({d % dv})'); continue
        step = -d // dv
        w = list(base); w[28730] += d; w[u] += step
        fwd(w, rounds=8)
        aw = L.all_atom_values(w)
        s = score(w)
        nzz = [a for a in range(L.NA) if aw[a]]
        print(f'  pair with x_{u}: a22230 {"0" if aw[22230]==0 else "nz"} '
              f'a22231 {"0" if aw[22231]==0 else "nz"} a37887 '
              f'{"0" if aw[37887]==0 else "nz"} a7930 {"0" if aw[7930]==0 else "nz"}'
              f'  nonzero {nzz}  score {s}')
        if s > 39026:
            T.save(w, os.path.join(HERE, f'D2_{s}.json'))
            print(f'    *** BEATS THE DELIVERABLE -- saved D2_{s}.json')
