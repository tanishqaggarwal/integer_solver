#!/usr/bin/env python3
"""Solve variables = free-input handles + the 220 wire members (as columns). Base at wire=1.
Build the full global Jacobian, test consistency of J*delta = -residual mod p (drive all 39033
roots to 0), and if consistent produce delta."""
import json, pickle, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire)
env.forced = {v: (s % p) for v, s in wire.items()}   # wire=1
env.jac_free = wireset                                 # treat wire members as solve columns
env.set_from_solution(best)
env.tangent_linear()
res = env.all_root_residuals()
print(f"[e] base wire=1: failing={len(res)}; solve columns = {len(env.freeinp)} handles + {len(wireset)} wire")

rows = []
active = set(); nnz = 0
for i in range(len(env.root_poly)):
    g = env.root_grad(i)
    if g:
        rows.append((i, g))
        for c in g: active.add(c)
        nnz += len(g)
wire_active = active & wireset
print(f"[e] nonzero rows={len(rows)}, active cols={len(active)} (wire active={len(wire_active)}), nnz={nnz}")

def rref(rows, res, track=True):
    pivots = {}; coldeg = defaultdict(int)
    for _, rd in rows:
        for c in rd: coldeg[c] += 1
    incons = 0; ilist = []; piv_order = []
    for k in sorted(range(len(rows)), key=lambda k: len(rows[k][1])):
        i, rd0 = rows[k]; rd = dict(rd0); rhs = (-res[i]) % p if i in res else 0
        while True:
            pc = None
            for c in rd:
                if c in pivots: pc = c; break
            if pc is None: break
            f = rd[pc]; prow, prhs = pivots[pc]
            for c, v in prow.items():
                nv = (rd.get(c,0)-f*v) % p
                if nv: rd[c] = nv
                elif c in rd: del rd[c]
            rhs = (rhs - f*prhs) % p
        if not rd:
            if rhs % p: incons += 1; ilist.append(i)
            continue
        pc = min(rd, key=lambda c: coldeg.get(c,0)); inv = pow(rd[pc],p-2,p)
        pivots[pc] = ({c:(v*inv)%p for c,v in rd.items()}, (rhs*inv)%p); piv_order.append(pc)
    delta = {}
    if track:
        for pc in reversed(piv_order):
            prow, prhs = pivots[pc]; s = prhs
            for c, v in prow.items():
                if c != pc:
                    dv = delta.get(c,0)
                    if dv: s = (s - v*dv) % p
            if s: delta[pc] = s
    return pivots, len(pivots), incons, ilist, delta

t0 = time.time()
pivots, rank, incons, ilist, delta = rref(rows, res)
print(f"[e] rank={rank}, INCONSISTENT rows={incons} {ilist[:10]}  ({time.time()-t0:.1f}s)")
print(f"[e] => (handles+wire) linear system is {'CONSISTENT' if incons==0 else 'INCONSISTENT'}")

if incons == 0:
    # apply delta: handles are free inputs, wire members via forced
    print(f"[e] applying delta over {len(delta)} columns...")
    for c, d in delta.items():
        if c in wireset:
            env.forced[c] = (env.forced[c] + d) % p
        else:
            env.valp[c] = (env.valp[c] + d) % p
    env.forward()
    res2 = env.all_root_residuals()
    print(f"[e] after 1 linear step: failing mod p = {len(res2)}  {sorted(res2)[:20]}")
    # save wire values reached
    pickle.dump({'delta': delta, 'wire_vals': {v: env.forced[v] for v in wireset}},
                open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire1e.pkl','wb'))
    print("[e] saved delta + wire values")
