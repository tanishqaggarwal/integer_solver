"""S10 step 108: build the frame-2 repair BY HAND and measure it.

All seven residual checks are zeroable by choice of the five detached parameters
and the solo handles:
    p | x_9118   -> x_29854 = p*x_1329 = 5113045*x_9118
    p | x_8731   -> x_31864 = p*x_10903 = -x_8731
    x_642  = p*x_17325 ;  x_7068 = x_2099 + 7376877*x_642
    x_28730 = p*x_9413
The only cost is what the consumers of x_9118, x_8731, x_7068, x_28730 do.
"""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import DETACH, definer, ORDER, FREE, CHECKS, fwd, score
P = ad.P
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
print(f'delivered witness in frame 2: score {score(base)}')

def build(v, k1off=0, k2off=0):
    v = list(v)
    k1 = v[9118] // P + k1off
    v[9118] = k1 * P
    v[1329] = 5113045 * v[7075] * k1
    k2 = v[8731] // P + k2off
    v[8731] = k2 * P
    v[31864] = -v[7075] * v[8731]
    v[10903] = v[31864] // P
    v[29854] = 5113045 * v[7075] * v[9118]
    v[642] = v[17325] * P
    v[28730] = v[9413] * P
    fwd(v, rounds=6)
    v[7068] = v[2099] + 7376877 * v[642]
    fwd(v, rounds=6)
    v[7068] = v[2099] + 7376877 * v[642]
    return v

v = build(base)
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
fail = L.failing_eqs(av)
print(f'\nconstructed: score {L.NEQ - len(fail)}  failing {len(fail)}')
print(f'  nonzero atoms ({len(nz)}): {nz}')
for a in [22229, 22230, 35758, 35759, 35760, 35761, 35762]:
    print(f'    a{a}: {"ZERO" if av[a] == 0 else "nonzero"}')
print(f'  failing equations: {sorted(fail)[:30]}')
for a in nz:
    print(f'    a{a:<6} ({len(L.atom2eq[a]):>2} eqs) '
          f'{"CHECK" if a in CHECKS else "gate"}')
if len(fail) < 7:
    T.save(v, os.path.join(HERE, f'construct_{L.NEQ-len(fail)}.json'))
    print(f'  saved construct_{L.NEQ-len(fail)}.json')
