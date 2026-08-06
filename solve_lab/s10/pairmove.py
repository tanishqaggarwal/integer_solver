"""S11 step 30: the compensating pair  x_28730 += d , x_19964 -= d.

a22231 = x_4432 - x_19964 - x_28730 and a37887 depends on x_4432 and on
(x_19964 + x_28730) only, so this move leaves BOTH untouched -- and x_4432 keeps
its value, so a7930 never sees it.  Meanwhile a22230 = x_28730 - p*x_9413 moves
freely, so A1 becomes free and eq 29125 is bought.  Cost = x_19964's other
consumers.
"""
import os, sys, random
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import definer, ORDER, FREE, CHECKS, fwd, score
P = ad.P
SSET = {22229, 22230, 35758, 35759, 35760, 35761, 35762}
random.seed(17)
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
print(f'frame 2 base {score(base)}')
print(f'x_19964 free? {19964 in FREE}; atoms {sorted(L.var_atoms[19964])}')
print(f'x_28730 atoms {sorted(L.var_atoms[28730])}')
print(f'coefficient of x_4432*x_19964 in a37887: '
      f'{L.polys[37887].get(tuple(sorted((4432,19964))), 0)}')
print(f'coefficient of x_4432*x_28730 in a37887: '
      f'{L.polys[37887].get(tuple(sorted((4432,28730))), 0)}')

print('\ncost of the compensating pair (atoms broken outside the seven):')
for lbl, d in (('+1', 1), ('+p', P), ('+rand', random.randrange(1, P))):
    v = list(base)
    v[28730] = v[28730] + d
    v[19964] = v[19964] - d
    fwd(v, rounds=8)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a] and a not in SSET]
    eqs = set()
    for a in nz: eqs |= set(L.atom2eq[a])
    print(f'  d = {lbl:<8}: outside-seven {nz} ({len(eqs)} eqs)  '
          f'a22231 {"0" if av[22231]==0 else "nz"}  a37887 '
          f'{"0" if av[37887]==0 else "nz"}  score {score(v)}')

print('\n=== now set A1 = 0 with the pair ===')
v = list(base)
d = v[9413] * P - v[28730]                # so that x_28730 becomes p*x_9413
v[28730] = v[28730] + d
v[19964] = v[19964] - d
fwd(v, rounds=8)
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
print(f'a22230 = {av[22230]}   a22231 = {"0" if av[22231]==0 else "nz"}   '
      f'a37887 = {"0" if av[37887]==0 else "nz"}')
print(f'nonzero {nz}  failing {len(L.failing_eqs(av))}  score {score(v)}')
if score(v) > 39026:
    T.save(v, os.path.join(HERE, f'PAIR_{score(v)}.json'))
    print(f'  *** BEATS THE DELIVERABLE -- saved PAIR_{score(v)}.json')
