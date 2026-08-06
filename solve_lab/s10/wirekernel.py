"""S10 step 50: THE decisive computation.

Write each wire member as  w_u = p + d_u.  Then every wire-identity atom is
linear and homogeneous in d:

    copy atom  x_i - x_j   ->  d_i - d_j
    root pin   x_26064 - p ->  d_root

Every equation containing such an atom becomes a homogeneous linear form in
d in Z^220.  Let M be that matrix.

    d = 0            -> the original p-wire (all equations satisfied, handles p-quantised)
    d in ker(M), d!=0 -> a wire deviation that costs NOTHING and unquantises handles

So: compute ker(M).  If some kernel vector has a coordinate that can be driven to
1 - p (making that wire member 1), the trapdoor's p-quantisation collapses for free.
"""
import os, sys, collections, json
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
print(f'wire members: {N}')

# --- identify the wire-identity atoms and their linear form in d -------------
IDENT = {}          # atom -> dict(wire index -> coeff)
for a in range(L.NA):
    vs = L.avars[a]
    if not vs or not all(u in widx for u in vs):
        continue
    poly = L.polys[a]
    form = collections.defaultdict(int)
    ok = True
    const = 0
    for m, c in poly.items():
        if len(m) == 0:
            const += c
        elif len(m) == 1:
            form[widx[m[0]]] += c
        else:
            ok = False; break
    if not ok:
        continue
    IDENT[a] = dict(form)
print(f'wire-identity atoms (linear, all variables in the wire): {len(IDENT)}')
print(f'  includes root pin a37694: {37694 in IDENT}')

# --- equations containing them ---------------------------------------------
EQS = set()
for a in IDENT:
    EQS |= set(L.atom2eq.get(a, ()))
EQS = sorted(EQS)
print(f'equations containing a wire-identity atom: {len(EQS)}')

rows = []
for e in EQS:
    m, sq, co = L.eq_atoms[e]
    form = collections.defaultdict(int)
    for a, c in co.items():
        if a in IDENT:
            for j, cc in IDENT[a].items():
                form[j] += c * cc
    form = {j: c for j, c in form.items() if c}
    if form:
        rows.append((e, form))
print(f'non-trivial rows: {len(rows)}')

# --- rank / kernel over Q ---------------------------------------------------
M = [[r.get(j, 0) for j in range(N)] for _, r in rows]
m = len(M)
A = [[Fraction(x) for x in row] for row in M]
piv = []
r = 0
for c in range(N):
    k = next((i for i in range(r, m) if A[i][c] != 0), None)
    if k is None:
        continue
    A[r], A[k] = A[k], A[r]
    pv = A[r][c]
    A[r] = [x / pv for x in A[r]]
    for i in range(m):
        if i != r and A[i][c] != 0:
            f = A[i][c]
            A[i] = [A[i][j] - f * A[r][j] for j in range(N)]
    piv.append(c); r += 1
    if r == m:
        break
rank = r
print(f'\nrank(M) = {rank} of {N} columns')
print(f'KERNEL DIMENSION = {N - rank}')

free_cols = [c for c in range(N) if c not in piv]
if not free_cols:
    print('\n*** kernel is trivial: the wire is RIGID -- w = p exactly, no free deviation')
else:
    print(f'\n*** NONTRIVIAL KERNEL: {len(free_cols)} free directions')
    basis = []
    for fc in free_cols:
        vec = [Fraction(0)] * N
        vec[fc] = Fraction(1)
        for i, c in enumerate(piv):
            vec[c] = -A[i][fc]
        den = 1
        for x in vec:
            den = den * x.denominator // __import__('math').gcd(den, x.denominator)
        ivec = [int(x * den) for x in vec]
        g = 0
        for x in ivec:
            g = __import__('math').gcd(g, x)
        if g:
            ivec = [x // g for x in ivec]
        basis.append(ivec)
    print(f'   kernel basis vectors: {len(basis)}')
    for b in basis[:5]:
        nzc = [(WIRE[j], b[j]) for j in range(N) if b[j]]
        print(f'   support {len(nzc)}: {nzc[:12]}')
    json.dump({'wire': WIRE, 'basis': basis},
              open(os.path.join(HERE, 'wirekernel.json'), 'w'))
    print('   saved wirekernel.json')
