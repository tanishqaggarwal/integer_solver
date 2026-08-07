"""S11 step 16: how expensive is it to MOVE rho = R1/R2 at all?

rho depends only on x_9118 and x_8731 mod p.  Measure, in frame 2, the exact cost
of perturbing each: which checks break, and are they cheap?
"""
import os, sys, random
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import definer, ORDER, FREE, CHECKS, fwd, score, grad
P = ad.P
random.seed(3)
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
av0 = L.all_atom_values(base)
SEVEN = {22229, 22230, 35758, 35759, 35760, 35761, 35762}
print(f'delivered witness {score(base)}; nonzero {[a for a in range(L.NA) if av0[a]]}')
vm = [x % P for x in base]
for a in (21617, 29539, 7930, 33796):
    g = grad(a, vm)
    print(f'  a{a}: d/dx_9118 = {"NONZERO" if g.get(9118,0)%P else "zero"}, '
          f'd/dx_8731 = {"NONZERO" if g.get(8731,0)%P else "zero"}  '
          f'(support {len(g)})')

print('\ncost of perturbing the two rho knobs (frame 2, seven checks ignored):')
for u in (9118, 8731):
    for lbl, d in (('+1', 1), ('+random', random.randrange(1, P))):
        v = list(base); v[u] = v[u] + d
        fwd(v, rounds=8)
        av = L.all_atom_values(v)
        nz = [a for a in range(L.NA) if av[a] and a not in SEVEN]
        eqs = set()
        for a in nz: eqs |= set(L.atom2eq[a])
        print(f'  x_{u} {lbl:<8}: broken non-seven atoms {len(nz)} {nz[:10]}')
        print(f'      equations they touch: {len(eqs)}')
