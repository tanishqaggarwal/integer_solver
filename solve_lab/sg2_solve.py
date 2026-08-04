#!/usr/bin/env python3
"""Solve the mod-p linear system: keep 602 checks=0, fix G1,G2 residues.
Sparse Gaussian elimination mod p (RREF). Report consistency and a particular delta."""
import pickle, time
p = 2**256-2**32-977
J = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/jac.pkl','rb'))
J_rows = J['J_rows']; constraints = J['constraints']; res = J['res']; ci = J['ci']
G1, G2 = 20862, 20864
def inv(a): return pow(a % p, p-2, p)

# Each row: (d=dict{col:coef}, b=rhs, ai)
rows = []
for ai in constraints:
    d = {c:v%p for c,v in J_rows[ci[ai]].items() if v%p}
    b = ((-res[ai]) % p) if ai in (G1,G2) else 0
    rows.append((d, b, ai))
print(f"rows: {len(rows)}")

# pivots: col -> (rowdict without pivot col (coef1 implied), rhs)
pivots = {}
def reduce_row(d, b):
    d = dict(d)
    # repeatedly eliminate any pivot column present
    stack = [c for c in d if c in pivots]
    while stack:
        col = stack.pop()
        if col not in d: continue
        f = d.pop(col)  # pivot coef is 1 (implied); remove it
        pd, pb = pivots[col]
        for c2, v2 in pd.items():
            nv = (d.get(c2,0) - f*v2) % p
            if nv:
                if c2 not in d and c2 in pivots: stack.append(c2)
                d[c2] = nv
            else:
                d.pop(c2, None)
        b = (b - f*pb) % p
    return d, b

t0 = time.time()
inconsistent = []
for d, b, ai in rows:
    d, b = reduce_row(d, b)
    if not d:
        if b % p != 0: inconsistent.append((ai, b))
        continue
    pcol = min(d.keys())
    ip = inv(d[pcol])
    nd = {c:(v*ip)%p for c,v in d.items() if c != pcol}
    nb = (b*ip)%p
    # back-reduce existing pivots containing pcol
    for c0 in list(pivots):
        pd, pb = pivots[c0]
        if pcol in pd:
            f = pd.pop(pcol)
            for c2,v2 in nd.items():
                nv = (pd.get(c2,0) - f*v2)%p
                if nv: pd[c2]=nv
                else: pd.pop(c2,None)
            pivots[c0] = (pd, (pb - f*nb)%p)
    pivots[pcol] = (nd, nb)
print(f"elimination done in {time.time()-t0:.1f}s  rank={len(pivots)}")
print(f"INCONSISTENT rows: {len(inconsistent)}")
for ai,b in inconsistent[:10]:
    print(f"   atom {ai}: 0 = {b}")
if not inconsistent:
    print(">>> CONSISTENT mod p. Missing DOF EXISTS.")
    delta = {pcol: pb for pcol,(pd,pb) in pivots.items()}  # free cols=0
    nz = {c:v for c,v in delta.items() if v%p}
    print(f"particular delta: {len(nz)} nonzero free-input residue changes")
    pickle.dump({'delta_modp':delta,'pivots':pivots}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/delta.pkl','wb'))
    print("saved delta.pkl")
