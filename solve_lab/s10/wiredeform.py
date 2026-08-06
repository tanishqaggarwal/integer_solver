"""S10 step 51: exact integer kernel of the wire deformation, and what values each
wire member can take for free.

ker has dimension 3 over Q.  Compute a saturated Z-basis; then the achievable set
of d_u is g_u * Z with g_u = gcd of the u-th coordinates over that basis, so

        w_u  in  p + g_u * Z.

If g_u divides p - 1 (or p + 1) we can set w_u = 1 (or -1) FOR FREE, and every
handle multiplied by that wire member becomes unquantised.
"""
import os, sys, collections, json, math
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
        rows.append([form.get(j, 0) for j in range(N)])
print(f'{len(rows)} rows x {N} cols')


def int_kernel(mat):
    m, n = len(mat), len(mat[0])
    A = [r[:] for r in mat]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    piv = []
    for r in range(m):
        while True:
            nz = [c for c in range(n) if c not in piv and A[r][c] != 0]
            if len(nz) <= 1: break
            nz.sort(key=lambda c: abs(A[r][c])); p0 = nz[0]
            for c in nz[1:]:
                q = A[r][c] // A[r][p0]
                if q:
                    for i in range(m): A[i][c] -= q * A[i][p0]
                    for i in range(n): U[i][c] -= q * U[i][p0]
        nz = [c for c in range(n) if c not in piv and A[r][c] != 0]
        if nz: piv.append(nz[0])
    return [[U[i][c] for i in range(n)] for c in range(n) if c not in piv]


B = int_kernel(rows)
print(f'integer kernel rank: {len(B)}')
# verify
for b in B:
    for r in rows:
        assert sum(r[j] * b[j] for j in range(N)) == 0
print('verified: every basis vector annihilates every row')

pm1, pp1 = P - 1, P + 1
hits = []
for j in range(N):
    g = 0
    for b in B:
        g = math.gcd(g, abs(b[j]))
    if g == 0:
        continue
    rec = {'var': WIRE[j], 'g_digits': len(str(g))}
    if pm1 % g == 0:
        rec['w=1'] = True
    if pp1 % g == 0:
        rec['w=-1'] = True
    if g == 1:
        rec['ANY'] = True
    if any(k in rec for k in ('w=1', 'w=-1', 'ANY')):
        hits.append(rec)
print(f'\nwire members with a nonzero free direction: '
      f'{sum(1 for j in range(N) if any(b[j] for b in B))} of {N}')
print(f'members where w = +-1 is reachable for free: {len(hits)}')
for h in hits[:20]:
    print('   ', h)

# smallest reachable |w| overall
best = []
for j in range(N):
    g = 0
    for b in B:
        g = math.gcd(g, abs(b[j]))
    if g == 0:
        continue
    r = P % g
    best.append((min(r, g - r), WIRE[j], g))
best.sort()
print('\nsmallest reachable |w_u| (want 1):')
for sm, u, g in best[:12]:
    print(f'   x_{u:<7} min|w| has {len(str(sm))} digits, g has {len(str(g))} digits')
json.dump({'wire': WIRE, 'basis': B}, open(os.path.join(HERE, 'wirebasis.json'), 'w'))
