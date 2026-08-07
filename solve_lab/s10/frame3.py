"""S11 step 27: frame 3 = frame 2 + detach x_4432 (a22231 becomes a check).

Then x_28730 no longer drives x_4432, so it stops feeding a7930; and K = x_4432 -
x_19964 (mod p) becomes tunable, which is exactly the compatibility condition the
1-dimensional kernel needs.  Price both knobs.
"""
import os, sys, random, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
DETACH = {7068: 22229, 28730: 22230, 29854: 35758, 31864: 35761, 642: 35762,
          4432: 22231}
definer = {t: a for t, a in L.definer.items() if t not in DETACH}
atom_out = {a: o for a, o in L.atom_out.items() if a not in set(DETACH.values())}
ORDER = [t for t in ad.ORDER if t not in DETACH]
FREE = set(t for t in range(L.NVARS) if t not in definer)
CHECKS = [a for a in range(L.NA) if a not in atom_out]
SSET = {22229, 22230, 35758, 35759, 35760, 35761, 35762, 22231}
random.seed(5)

def fwd(v, rounds=6):
    for _ in range(rounds):
        for u in ORDER:
            nv = T.solve_lin(definer[u], u, v)
            if nv is not None: v[u] = nv
    return v
def score(v):
    return L.NEQ - len(L.failing_eqs(L.all_atom_values(v)))

if __name__ == '__main__':
    base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
    b2 = list(base); fwd(b2)
    print(f'delivered witness in frame 3: {score(b2)} (on-manifold: {score(b2)==39026})')
    av = L.all_atom_values(b2)
    print(f'  nonzero: {[a for a in range(L.NA) if av[a]]}')
    print(f'  a22231 = {av[22231]}')
    print('\nprice of each knob (atoms broken OUTSIDE the eight):')
    for u in (4432, 28730, 7068, 9118, 8731):
        for lbl, d in (('+1', 1), ('+rand', random.randrange(1, P))):
            v = list(b2); v[u] = v[u] + d
            fwd(v, rounds=8)
            a2 = L.all_atom_values(v)
            nz = [a for a in range(L.NA) if a2[a] and a not in SSET]
            eqs = set()
            for a in nz: eqs |= set(L.atom2eq[a])
            print(f'  x_{u:<6} {lbl:<6}: outside-eight {nz}  ({len(eqs)} equations)')
