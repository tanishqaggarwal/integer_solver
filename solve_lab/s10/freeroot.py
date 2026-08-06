"""S10 step 74: construct the deformation that frees the wire root for 1 equation.

rootfree.py: e_root = y0^T M with supp(y0) = {eq 37257}.  So there is a wire
deformation satisfying ALL identity equations except 37257, with d_root != 0 --
the root pin a37694 is compensated by copy atoms in 11 of its 12 equations, and
only 37257 has no other identity atom to absorb it.

Compute ker(M without row 37257) (dim 4), extract a direction with d_root != 0,
apply it, and measure.  If the square-check damage is small the trapdoor falls.
"""
import os, sys, json, collections, math
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
base = L.load(os.path.join(HERE, 'forward_state.json'))
WIRE = sorted(u for u in range(L.NVARS) if base[u] == P)
widx = {u: i for i, u in enumerate(WIRE)}
N = len(WIRE)
ROOT = widx[26064]
DROP = 37257

IDENT = {}
for a in range(L.NA):
    vs = L.avars[a]
    if not vs or not all(u in widx for u in vs):
        continue
    form = collections.defaultdict(int); ok = True
    for m, c in L.polys[a].items():
        if len(m) == 1:
            form[widx[m[0]]] += c
        elif len(m) != 0:
            ok = False; break
    if ok:
        IDENT[a] = dict(form)
EQS = sorted(set().union(*[set(L.atom2eq.get(a, ())) for a in IDENT]))
rows, rowid = [], []
for e in EQS:
    m, sq, co = L.eq_atoms[e]
    form = collections.defaultdict(int)
    for a, c in co.items():
        if a in IDENT:
            for j, cc in IDENT[a].items():
                form[j] += c * cc
    form = {j: c for j, c in form.items() if c}
    if form:
        rows.append([form.get(j, 0) for j in range(N)]); rowid.append(e)

print(f'eq {DROP} identity-form nonzero coords: '
      f'{[j for j in range(N) if rows[rowid.index(DROP)][j]]}  '
      f'(root index {ROOT})')

keep = [rows[i] for i in range(len(rows)) if rowid[i] != DROP]
print(f'dropping eq {DROP}: {len(keep)} identity equations remain')


def kernel_q(mat, n):
    m = len(mat)
    A = [[Fraction(x) for x in r] for r in mat]
    piv = []; r = 0
    for c in range(n):
        k = next((i for i in range(r, m) if A[i][c] != 0), None)
        if k is None: continue
        A[r], A[k] = A[k], A[r]
        pv = A[r][c]; A[r] = [x / pv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(n)]
        piv.append(c); r += 1
        if r == m: break
    free = [c for c in range(n) if c not in piv]
    basis = []
    for fc in free:
        vec = [Fraction(0)] * n
        vec[fc] = Fraction(1)
        for i, c in enumerate(piv):
            vec[c] = -A[i][fc]
        den = 1
        for x in vec:
            den = den * x.denominator // math.gcd(den, x.denominator)
        iv = [int(x * den) for x in vec]
        g = 0
        for x in iv: g = math.gcd(g, x)
        if g: iv = [x // g for x in iv]
        basis.append(iv)
    return basis, r


B, rk = kernel_q(keep, N)
print(f'rank = {rk}; kernel dimension = {len(B)}')
withroot = [b for b in B if b[ROOT] != 0]
print(f'kernel directions with d_root != 0: {len(withroot)} of {len(B)}')
for i, b in enumerate(B):
    nz = [j for j in range(N) if b[j]]
    print(f'  dir {i}: support {len(nz)}  d_root = '
          f'{"0" if b[ROOT] == 0 else str(b[ROOT])[:30] + "... (" + str(len(str(abs(b[ROOT])))) + " digits)"}')

if not withroot:
    print('\nno direction moves the root -- contradicts rootfree.py, re-check')
    sys.exit()

WIREDEF = set(L.definer[u] for u in WIRE if u in L.definer)


def fwd(v, rounds=3):
    for _ in range(rounds):
        for x in ad.ORDER:
            a = L.definer[x]
            if a in WIREDEF:
                continue
            nv = T.solve_lin(a, x, v)
            if nv is not None:
                v[x] = nv
    return v


f0 = len(L.failing_eqs(L.all_atom_values(base)))
print(f'\nbase failing={f0} score={L.NEQ-f0}')
for i, b in enumerate(B):
    if b[ROOT] == 0:
        continue
    for s in (1, -1):
        v = list(base)
        for j, u in enumerate(WIRE):
            v[u] += s * b[j]
        fwd(v)
        av = L.all_atom_values(v)
        nz = [a for a in range(L.NA) if av[a]]
        hard = [a for a in nz if a not in WIREDEF]
        fail = L.failing_eqs(av)
        print(f'  dir {i} x{s:+d}: w_root now {len(str(abs(v[26064])))} digits; '
              f'nonzero atoms {len(nz)} (non-copy {len(hard)}); '
              f'failing {len(fail)} score {L.NEQ-len(fail)}', flush=True)
        print(f'      non-copy broken: {hard[:16]}')
        if len(fail) < f0:
            T.save(v, os.path.join(HERE, f'freeroot_{L.NEQ-len(fail)}.json'))
