"""S11 step 12: zero SIX of the seven, leave a35758 nonzero.

lattice11.py: with only a35758 nonzero exactly 6 equations fail -> 39,027.
That configuration needs a35759 = 0, i.e. x_29854 = 5113045*x_7075*x_9118 with
x_29854 detached -- so x_9118 is NOT forced to be a multiple of p.  My earlier
construct2.py imposed BOTH p|x_9118 and p|x_8731; letting a35758 break is exactly
what buys back x_9118.  Start from the delivered witness, where the cluster is
already satisfied.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import DETACH, definer, ORDER, FREE, CHECKS, fwd, score
P = ad.P
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
print(f'delivered witness in frame 2: {score(base)}')

def build(v, k2off=0, x1329=None):
    v = list(v)
    # a35760 = a35761 = 0  requires p | x_8731
    k2 = v[8731] // P + k2off
    v[8731] = k2 * P
    v[31864] = -v[7075] * v[8731]
    v[10903] = v[31864] // P
    # a35759 = 0 through the DETACHED x_29854 ; x_9118 left alone
    v[29854] = 5113045 * v[7075] * v[9118]
    # a35758 = x_29854 - p*x_1329 : deliberately NONZERO
    v[1329] = v[29854] // P if x1329 is None else x1329
    v[642] = v[17325] * P
    v[28730] = v[9413] * P
    fwd(v, rounds=6)
    v[7068] = v[2099] + 7376877 * v[642]
    fwd(v, rounds=6)
    v[7068] = v[2099] + 7376877 * v[642]
    return v

for k2 in (0, 1, -1):
    v = build(base, k2off=k2)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    print(f'\nk2off={k2}: failing {len(fail)}  score {L.NEQ-len(fail)}  nonzero {nz}')
    for a in (22229, 22230, 35758, 35759, 35760, 35761, 35762):
        print(f'   a{a}: {"ZERO" if av[a]==0 else "NONZERO"}', end='')
    print()
    if len(fail) < 7:
        T.save(v, os.path.join(HERE, f'C3_{L.NEQ-len(fail)}.json'))
        print(f'   *** saved C3_{L.NEQ-len(fail)}.json')
