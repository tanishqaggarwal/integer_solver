"""S10 step 45: joint integer solve of the four coupled checks on the wire=1 branch.

a7930's handle x_11052 also sits in a41512; a29539's handle x_30163 also sits in
a40826.  All four are broken, so solve them together as one integer linear system
in the free inputs / solo handles they share.
"""
import os, sys, collections, math, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
ROOTATOM = 37694
BLOCK = {ROOTATOM}
TARGET = [7930, 29539, 39417, 40826, 41512]

v = L.load(os.path.join(HERE, 'wire1_solved2.json'))
av = L.all_atom_values(v)
print('current nonzero:', [a for a in range(L.NA) if av[a]])
for a in TARGET:
    print(f'  a{a} = {str(av[a])[:44]}  (neq={len(L.atom2eq.get(a,{}))})')


def fwd_block(v, rounds=4):
    for _ in range(rounds):
        for u in ad.ORDER:
            d = L.definer[u]
            if d in BLOCK:
                continue
            nv = T.solve_lin(d, u, v)
            if nv is not None:
                v[u] = nv
    return v


# candidate knobs: free inputs inside the target atoms, plus solo handles of the
# defining atoms of any gate output inside them
cand = set()
for a in TARGET:
    for u in L.avars[a]:
        if u not in L.definer:
            cand.add(u)
        else:
            d = L.definer[u]
            for h in L.avars[d]:
                if h not in L.definer:
                    cand.add(h)
cand = sorted(cand)
print(f'\ncandidate knobs: {len(cand)}')

# exact effect of each knob on each target atom (finite difference of +1, exact
# because everything here is affine in these variables once forward-evaluated)
base = {a: av[a] for a in TARGET}
J = {}
safe = []
for u in cand:
    w = list(v); w[u] = w[u] + 1
    fwd_block(w)
    aw = L.all_atom_values(w)
    col = {a: aw[a] - base[a] for a in TARGET}
    extra = sorted(set(x for x in range(L.NA) if aw[x]) - set(x for x in range(L.NA) if av[x]))
    if extra:
        continue                       # knob has collateral; skip
    if any(col.values()):
        J[u] = col
        safe.append(u)
print(f'collateral-free knobs that move a target: {len(safe)} -> {[f"x_{u}" for u in safe]}')
for u in safe:
    print(f'   x_{u:<7} {[(a, str(d)[:22]) for a, d in J[u].items() if d]}')

# integer linear solve  sum_u y_u * J[u][a]  =  -base[a]   for all a in TARGET
rows = [[J[u].get(a, 0) for u in safe] for a in TARGET]
rhs = [-base[a] for a in TARGET]


def solve_int(rows, rhs, n):
    m = len(rows)
    A = [r[:] for r in rows]
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
        piv.append(nz[0] if nz else None)
    w = [0] * n; b = list(rhs)
    for r in range(m):
        c = piv[r]
        if c is None:
            if b[r] != 0: return None
            continue
        if b[r] % A[r][c]: return None
        w[c] = b[r] // A[r][c]
        for rr in range(r + 1, m): b[rr] -= A[rr][c] * w[c]
    return [sum(U[i][c] * w[c] for c in range(n)) for i in range(n)]


sol = solve_int(rows, rhs, len(safe)) if safe else None
if sol is None:
    print('\n*** joint integer system has no solution over these knobs')
else:
    print(f'\n*** solution found: {[(f"x_{safe[i]}", str(sol[i])[:24]) for i in range(len(safe)) if sol[i]]}')
    for i, u in enumerate(safe):
        v[u] = v[u] + sol[i]
    fwd_block(v)
    av2 = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av2[a]]
    fail = L.failing_eqs(av2)
    print(f'   after: nz={nz} failing={len(fail)} score={L.NEQ-len(fail)}')
    T.save(v, os.path.join(HERE, 'joint_out.json'))
