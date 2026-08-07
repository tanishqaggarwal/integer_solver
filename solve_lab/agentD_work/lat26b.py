"""Exact max-satisfiable over the 12 equations of the 39,026 placement,
using the 6 cost-free generators measured with the broken gates BLOCKED."""
import json, itertools, sys
import dlib as L
from intsolve import solve_int

d = json.load(open('gens26.json'))
COLS = d['cols']
base = [int(x) for x in d['base']]
gens = [(u, [int(x) for x in g]) for u, g in d['gens']]
print('cols', COLS)
print('gens', [u for u, _ in gens])

Eq = set()
for a in (22229, 22230, 35758, 35759, 35760, 35761, 35762):
    Eq |= set(L.atom2eq.get(a, {}))
Eq = sorted(Eq)
M = []
for i in Eq:
    m, sq, co = L.eq_atoms[i]
    M.append([co.get(a, 0) for a in COLS])
print('E', len(Eq), Eq)

G = [[g[j] for u, g in gens] for j in range(len(COLS))]
ng = len(gens)

for size in range(len(Eq), 0, -1):
    found = None
    cnt = 0
    for S in itertools.combinations(range(len(Eq)), size):
        A = [[sum(M[si][j] * G[j][k] for j in range(len(COLS))) for k in range(ng)] for si in S]
        b = [-sum(M[si][j] * base[j] for j in range(len(COLS))) for si in S]
        x = solve_int(A, b)
        cnt += 1
        if x is not None:
            found = (S, x)
            break
    print(f'size {size}: tried {cnt} subsets, solvable={found is not None}', flush=True)
    if found:
        print('  S =', [Eq[i] for i in found[0]])
        print('  k =', found[1])
        print('MAX =', size, ' failing =', len(Eq) - size, ' score =', 39033 - (len(Eq) - size))
        json.dump({'S': [Eq[i] for i in found[0]], 'k': [str(z) for z in found[1]],
                   'gens': [u for u, _ in gens]}, open('lat26b.json', 'w'))
        break
