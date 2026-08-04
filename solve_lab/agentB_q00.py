#!/usr/bin/env python3
"""Quadrant (0,0): override={} (no activators), all free inputs = 0. Build Jacobian mod p at
all-zero, report failing roots, rank, and consistency of J*delta = -residual."""
import time, sys
from collections import defaultdict
from agentB_setup import build, load, Env, p, NVARS

CQ = '/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_cache_q00.pkl'
t0 = time.time()
import os
data = load(ov={}, cache=CQ) if os.path.exists(CQ) else build(ov={}, cache=CQ)
env = Env(data)
print(f"[q00] free inputs={len(env.freeinp)}  (override empty)", file=sys.stderr)
env.set_free({})   # all zero
res = env.all_root_residuals()
print(f"[q00] nonzero roots at all-zero: {len(res)}  in {time.time()-t0:.1f}s")
print("[q00] indices:", sorted(res)[:40])

env.tangent_linear()
# build rows
sat_rows = []; core_rows = []
nnz = 0; cols = set()
for i in range(len(env.root_poly)):
    g = env.root_grad(i)
    if not g: continue
    nnz += len(g); cols |= set(g)
    if i in res: core_rows.append((i, g, (-res[i]) % p))
    else: sat_rows.append((i, g, 0))
print(f"[q00] nonzero rows: sat={len(sat_rows)} core={len(core_rows)}; nnz={nnz}; active cols={len(cols)}")

def rref_track(rows):
    pivots = {}; coldeg = defaultdict(int)
    for rd, _ in rows:
        for c in rd: coldeg[c] += 1
    incons = 0; npiv = 0
    order = sorted(range(len(rows)), key=lambda k: len(rows[k][0]))
    for k in order:
        rd, rhs = dict(rows[k][0]), rows[k][1]
        while True:
            pcol = None
            for c in rd:
                if c in pivots: pcol = c; break
            if pcol is None: break
            f = rd[pcol]; prow, prhs = pivots[pcol]
            for c, v in prow.items():
                nv = (rd.get(c, 0) - f*v) % p
                if nv: rd[c] = nv
                elif c in rd: del rd[c]
            rhs = (rhs - f*prhs) % p
        if not rd:
            if rhs % p: incons += 1
            continue
        pc = min(rd, key=lambda c: coldeg.get(c, 0))
        inv = pow(rd[pc], p-2, p)
        pivots[pc] = ({c:(v*inv)%p for c,v in rd.items()}, (rhs*inv)%p)
        npiv += 1
    return pivots, npiv, incons

# rank of satisfied constraints
piv_sat, rank_sat, _ = rref_track([(g, 0) for _, g, _ in sat_rows])
print(f"[q00] rank(J_sat) = {rank_sat}")
# add core rows one at a time
piv = dict(piv_sat); indep = 0; incons = 0; absorbed = 0
def reduce_against(pivots, rd0, rhs0):
    rd = dict(rd0); rhs = rhs0
    while True:
        pcol = None
        for c in rd:
            if c in pivots: pcol = c; break
        if pcol is None: break
        f = rd[pcol]; prow, prhs = pivots[pcol]
        for c, v in prow.items():
            nv = (rd.get(c, 0) - f*v) % p
            if nv: rd[c] = nv
            elif c in rd: del rd[c]
        rhs = (rhs - f*prhs) % p
    return rd, rhs
for i, g, rhs in core_rows:
    rd, r = reduce_against(piv, g, rhs)
    if not rd:
        if r % p: incons += 1
        else: absorbed += 1
    else:
        pc = min(rd, key=lambda c: c); inv = pow(rd[pc], p-2, p)
        piv[pc] = ({c:(v*inv)%p for c,v in rd.items()}, (r*inv)%p); indep += 1
print(f"[q00] core rows: independent(satisfiable)={indep}, absorbed={absorbed}, INCONSISTENT={incons}")
print(f"[q00] {'FULLY CONSISTENT' if incons==0 else 'INCONSISTENT'} linear system (t={time.time()-t0:.1f}s)")
