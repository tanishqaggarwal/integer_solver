"""RESIDUE-AWARE PRICER.

Every handle in the file has the form p*(free variable), so a violated atom's value can be shifted
only by multiples of p: the achievable set is  v in v0 + p*Z^S,  where v0 is read off a REAL
construction.  Hence equation e can vanish only if  M_e . v0 = 0 (mod p)  -- a condition on the
residues, invisible to any function of the incidence matrix alone.  Among the equations passing
that test, a subset T can vanish simultaneously iff  M_T w = -(M_T v0)/p  has an integer solution.

price(cluster) = |E| - max|T|.
"""
import sys, os, json, collections, itertools
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from close4 import *
W = os.path.dirname(os.path.abspath(__file__)) + '/'
AE = {a: frozenset(L.atom2eq.get(a, {})) for a in range(L.NA)}
e2a = collections.defaultdict(set)
for a, es in AE.items():
    for e in es: e2a[e].add(a)
def inside(E):
    c = set()
    for e in E: c |= e2a[e]
    return sorted(a for a in c if AE[a] <= E)

def int_solvable(rows, rhs):
    """Does rows * w = rhs have an integer solution?  Hermite-style elimination over Z."""
    A = [list(r) + [b] for r, b in zip(rows, rhs)]
    n = len(A); m = len(rows[0]) if n else 0
    r = 0
    for c in range(m):
        piv = None
        for i in range(r, n):
            if A[i][c]: piv = i; break
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        # integer row reduction by repeated gcd steps
        changed = True
        while changed:
            changed = False
            for i in range(r + 1, n):
                if A[i][c]:
                    q = A[i][c] // A[r][c]
                    if q: A[i] = [x - q * y for x, y in zip(A[i], A[r])]
                    if A[i][c]:
                        A[r], A[i] = A[i], A[r]; changed = True
        for i in range(r):
            if A[i][c] and A[r][c] and A[i][c] % A[r][c] == 0:
                q = A[i][c] // A[r][c]
                A[i] = [x - q * y for x, y in zip(A[i], A[r])]
        r += 1
    for i in range(r, n):
        if all(x == 0 for x in A[i][:m]) and A[i][m] != 0: return False
    for i in range(r):
        row = A[i]
        nz = [x for x in row[:m] if x]
        if not nz:
            if row[m]: return False
            continue
        g = 0
        for x in nz:
            while x: g, x = x, g % x
        g = abs(g)
        if g and row[m] % g: return False
    return True

def price(E, S, v0, forced=None, cap=4096):
    """v0 = list of atom values (integers) at the construction, aligned with S."""
    E = sorted(E)
    M = [[L.eq_atoms[e][2].get(a, 0) for a in S] for e in E]
    res = [sum(c * x for c, x in zip(row, v0)) % P for row in M]
    cand = [i for i in range(len(E)) if res[i] == 0]
    blocked = len(E) - len(cand)
    best = 0; bestT = None
    idx = list(range(len(cand)))
    for k in range(len(cand), -1, -1):
        if best: break
        if len(list(itertools.combinations(idx, k))) > cap and k < len(cand): break
        for combo in itertools.combinations(idx, k):
            T = [cand[i] for i in combo]
            rows = [M[i] for i in T]
            rhs = [-sum(c * x for c, x in zip(M[i], v0)) // P for i in T]
            if not rows or int_solvable(rows, rhs):
                best = k; bestT = T; break
    return dict(nE=len(E), nS=len(S), residue_blocked=blocked,
                max_satisfiable=best, price=len(E) - best)

if __name__ == '__main__':
    print('=' * 72)
    print('CALIBRATION: the 39,026 deliverable cluster (true cost 7)')
    d = json.load(open(W + '../best/new_instance_partial_39026.json'))
    v = [0] * L.NVARS
    for k, val in d.items(): v[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
    av = L.all_atom_values(v)
    D = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
    E = frozenset().union(*[AE[a] for a in D]); S = inside(E)
    v0 = [av[a] for a in S]
    r = price(E, S, v0)
    print('   ', r)
    ok = (r['price'] == 7)
    print('    CALIBRATION %s  (structural relaxation gave 5; residue-aware gives %d)'
          % ('PASSES' if ok else 'FAILS', r['price']))
    print('    nonzero atoms of S at the deliverable:', [a for a, x in zip(S, v0) if x])
