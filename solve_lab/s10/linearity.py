"""S11 step 4: is the free-input -> check map LINEAR mod p at this stratum?

If it is, the 1655x707 closure is not a linearisation at all -- it is the exact
system, and its inconsistency is final for this stratum.  If it is not, Newton
from a different starting point can reach another root.
Test: compare the exact finite difference against grad*delta for LARGE random
delta, on many free inputs and many checks.
"""
import os, sys, random, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
random.seed(12345)
v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm0 = [x % P for x in v0]
av0 = L.all_atom_values(v0)
CHECKS = [a for a in range(L.NA) if a not in atom_out]
U = sorted(set(ad.grad(21617, vm0)) | set(ad.grad(29539, vm0)))
print(f'testing {len(U)} free inputs from the cluster gradient support')

lin = nonlin = 0
nonlin_u = collections.Counter()
for u in U[:40]:
    col = jac_column(u, v0, vm0, CHECKS)
    for trial in range(2):
        d = random.randrange(1, P)
        w = list(v0); w[u] = w[u] + d
        ad.fwd(w, rounds=8)
        aw = L.all_atom_values(w)
        for c in list(col)[:25]:
            pred = (av0[c] + col[c] * d) % P
            got = aw[c] % P
            if pred == got: lin += 1
            else:
                nonlin += 1; nonlin_u[u] += 1
print(f'\nlinear predictions correct: {lin}   WRONG: {nonlin}')
print(f'  free inputs showing nonlinearity: {len(nonlin_u)} '
      f'{dict(nonlin_u.most_common(8))}')

# and the decisive one: are the two cluster residues linear in the free inputs?
print('\ncluster residues under large random moves:')
for u in U[:12]:
    col = jac_column(u, v0, vm0, CHECKS)
    d = random.randrange(1, P)
    w = list(v0); w[u] = w[u] + d
    ad.fwd(w, rounds=8)
    aw = L.all_atom_values(w)
    ok = []
    for c in (21617, 29539):
        pred = (av0[c] + col.get(c, 0) * d) % P
        ok.append(pred == aw[c] % P)
    print(f'  x_{u:<7} a21617 linear {ok[0]}   a29539 linear {ok[1]}')
