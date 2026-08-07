"""Independent re-derivation of the optimality claim for the 39,026 witness.

1. E = equations touched by its 7 nonzero atoms;  the 12 x 7 coefficient matrix.
2. Measure the reachable lattice of the 7 atom values by scanning every free input.
3. For every subset S of E, decide whether some reachable atom vector zeroes S.
"""
import json, sys, time, itertools, collections
from fractions import Fraction
import dlib as L
import engine2 as E
P = L.P

st0 = E.St(L.load(sys.argv[1] if len(sys.argv) > 1 else '../best/new_instance_partial_39026.json'))
NZ = st0.nz()
print('score', st0.score, 'nz', NZ)
Eq = set()
for a in NZ:
    Eq |= set(L.atom2eq.get(a, {}))
Eq = sorted(Eq)
print('E =', len(Eq), Eq)
atoms = set()
for i in Eq:
    atoms |= set(L.eq_atoms[i][2])
print('atoms appearing in E:', len(atoms), sorted(atoms))
conf = [a for a in sorted(atoms) if set(L.atom2eq.get(a, {})) <= set(Eq)]
print('confined atoms (out=0):', len(conf), conf)

COLS = conf
M = []
for i in Eq:
    m, sq, co = L.eq_atoms[i]
    M.append([co.get(a, 0) for a in COLS])
print('matrix', len(M), 'x', len(COLS))


def rank(rows):
    rows = [[Fraction(x) for x in r] for r in rows]
    n = len(rows[0]) if rows else 0
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, len(rows)):
            if rows[i][c]:
                piv = i
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pv = rows[r][c]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f = rows[i][c] / pv
                for j in range(c, n):
                    rows[i][j] -= f * rows[r][j]
        r += 1
    return r


print('rank(full 12 x n) =', rank(M))

# --- measure reachable deltas of the confined atom values ---
base = [st0.av[a] for a in COLS]
gens = []
t0 = time.time()
for k, u in enumerate(sorted(L.freeset)):
    st = st0.clone()
    st.apply({u: st.v[u] + 1})
    # only accept moves that leave every atom OUTSIDE E unchanged
    ok = True
    for a in range(L.NA):
        if st.av[a] != st0.av[a] and a not in atoms:
            ok = False
            break
    if not ok:
        continue
    d = [st.av[a] - base[i] for i, a in enumerate(COLS)]
    if any(d):
        gens.append((u, d))
    if k % 2000 == 0:
        print('  scan', k, f'{time.time()-t0:.0f}s gens={len(gens)}', flush=True)
print('cost-free generators of the residual lattice:', len(gens))
for u, d in gens:
    print('   x_%-7d %s' % (u, [str(x)[:26] for x in d]))
json.dump({'E': Eq, 'cols': COLS, 'M': M,
           'base': [str(x) for x in base],
           'gens': [[u, [str(x) for x in d]] for u, d in gens]}, open('opt26.json', 'w'))
