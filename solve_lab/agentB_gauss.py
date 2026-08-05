#!/usr/bin/env python3
"""Sparse Gaussian elimination mod p of the full root-Jacobian.
Reports rank(J), rank([J|b]), consistency of J*delta = b  (b_i = -residual).
If consistent, computes a particular solution delta (a mod-p Newton step)."""
import pickle, time, sys
from collections import defaultdict
p = 2**256 - 2**32 - 977
SC = '/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_jac.pkl'

t0 = time.time()
D = pickle.load(open(SC, 'rb'))
rows_in = D['rows']; res = D['res']
# build augmented rows: (rowdict col->val, rhs). b = -residual
rows = []
for i, g in rows_in:
    rhs = (-res[i]) % p if i in res else 0
    rows.append([dict(g), rhs])
print(f"[gauss] {len(rows)} rows; core(nonzero rhs)={sum(1 for r in rows if r[1])}")

# column degrees (for pivot choice, reduce fill)
def solve(rows):
    pivots = {}          # col -> [rowdict(normalized pivot=1), rhs]
    coldeg = defaultdict(int)
    for rd, _ in rows:
        for c in rd: coldeg[c] += 1
    incons = 0
    incons_core = 0
    order = sorted(range(len(rows)), key=lambda k: len(rows[k][0]))
    npiv = 0
    for oi, k in enumerate(order):
        rd, rhs = rows[k]
        rd = dict(rd)
        # reduce against existing pivots
        while True:
            pcol = None
            for c in rd:
                if c in pivots: pcol = c; break
            if pcol is None: break
            f = rd[pcol]
            prow, prhs = pivots[pcol]
            for c, v in prow.items():
                nv = (rd.get(c, 0) - f * v) % p
                if nv: rd[c] = nv
                elif c in rd: del rd[c]
            rhs = (rhs - f * prhs) % p
        if not rd:
            if rhs % p != 0:
                incons += 1
                if res.get(list(res)[0]) is not None and (k < len(rows_in)):
                    pass
            continue
        # choose pivot col: min coldeg
        pc = min(rd, key=lambda c: coldeg.get(c, 0))
        inv = pow(rd[pc], p - 2, p)
        prow = {c: (v * inv) % p for c, v in rd.items()}
        prhs = (rhs * inv) % p
        pivots[pc] = [prow, prhs]
        npiv += 1
        if oi % 2000 == 0:
            print(f"  processed {oi}/{len(order)}  pivots={npiv}  incons={incons}  t={time.time()-t0:.1f}s", flush=True)
    return pivots, npiv, incons

pivots, rank, incons = solve(rows)
print(f"[gauss] rank(J) (pivots) = {rank}")
print(f"[gauss] inconsistent (0=nonzero) rows encountered = {incons}")
print(f"[gauss] {'CONSISTENT' if incons==0 else 'INCONSISTENT'}: J*delta=b  (t={time.time()-t0:.1f}s)")

# Determine which core equations cause inconsistency: solve WITHOUT core rhs to get rank(J),
# then check if adding each core row's rhs is in row space.
