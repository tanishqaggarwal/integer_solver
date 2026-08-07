"""S11 step 17: price the two CONSTRAINED atom values separately.

A2..A6 are free at zero cost (x_9118, x_8731, x_29854, x_31864, x_642 and their
solo handles).  Only A0 = a22229 (through x_7068) and A1 = a22230 (through
x_28730) carry collateral.  Rows:
    eq 29125 = A1          -> satisfied iff A1 = 0
    eq  2554 = A0 + 13*A1  -> satisfied iff A0 = -13*A1
With A0, A1 given, the other 10 equations are 10 linear conditions on the 5 free
values A2..A6, so 5 are satisfiable.  Hence
    A1 != 0, A0 != 0  -> 5 satisfied, 7 failing   (the delivered witness)
    A1 == 0           -> 6 satisfied, 6 failing   -> 39,027
    A0 == A1 == 0     -> all 12                   -> 39,033
So: what does it cost to set A1 = 0, i.e. x_28730 = p*x_9413 ?
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import definer, ORDER, FREE, CHECKS, fwd, score
P = ad.P
SEVEN = {22229, 22230, 35758, 35759, 35760, 35761, 35762}
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
av0 = L.all_atom_values(base)
print(f'delivered witness {score(base)}')
print(f'  A0 = a22229 = {str(av0[22229])[:30]}')
print(f'  A1 = a22230 = {str(av0[22230])[:30]}')
print(f'  x_28730 mod p = {str(av0[22230] % P)[:30]}')

print('\nconsumers of the two constrained variables:')
for u in (7068, 28730, 642, 29854, 31864):
    print(f'  x_{u}: atoms {sorted(L.var_atoms[u])}')

print('\ncost of moving each detached variable (outside the seven):')
for u in (7068, 28730, 642, 29854, 31864):
    v = list(base); v[u] = v[u] + 1
    fwd(v, rounds=8)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a] and a not in SEVEN]
    eqs = set()
    for a in nz: eqs |= set(L.atom2eq[a])
    print(f'  x_{u:<6} +1 : broken atoms {len(nz)} {nz[:8]}  equations {len(eqs)}')

print('\n=== set A1 = 0 : x_28730 = p * x_9413 ===')
v = list(base)
v[28730] = v[9413] * P
fwd(v, rounds=8)
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
print(f'  a22230 now {"ZERO" if av[22230]==0 else "nonzero"}')
print(f'  nonzero atoms {nz}')
print(f'  failing {len(L.failing_eqs(av))}  score {score(v)}')
