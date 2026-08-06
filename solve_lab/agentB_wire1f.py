#!/usr/bin/env python3
"""Iterated mod-p Newton over (handles + 220 wire members) from the wire=1 base. Consistent linear
system each step -> should converge (heals 13 unpacking + solves 20 core)."""
import json, pickle, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire)
env.forced = {v: (s % p) for v, s in wire.items()}
env.jac_free = wireset
env.set_from_solution(best)

def step():
    env.tangent_linear()
    res = env.all_root_residuals()
    rows = []
    for i in range(len(env.root_poly)):
        g = env.root_grad(i)
        if g: rows.append((i, g))
    coldeg = defaultdict(int)
    for _, rd in rows:
        for c in rd: coldeg[c] += 1
    pivots = {}; piv_order = []; incons = 0
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
            if rhs % p: incons += 1
            continue
        pc = min(rd, key=lambda c: coldeg.get(c,0)); inv = pow(rd[pc],p-2,p)
        pivots[pc] = ({c:(v*inv)%p for c,v in rd.items()}, (rhs*inv)%p); piv_order.append(pc)
    delta = {}
    for pc in reversed(piv_order):
        prow, prhs = pivots[pc]; s = prhs
        for c, v in prow.items():
            if c != pc:
                dv = delta.get(c,0)
                if dv: s = (s - v*dv) % p
        if s: delta[pc] = s
    return len(res), incons, delta

for it in range(25):
    t = time.time()
    nfail, incons, delta = step()
    for c, d in delta.items():
        if c in wireset: env.forced[c] = (env.forced[c] + d) % p
        else: env.valp[c] = (env.valp[c] + d) % p
    env.forward()
    res2 = env.all_root_residuals()
    print(f"  it {it}: fail {nfail}->{len(res2)}, incons={incons}, |delta|={len(delta)}, {time.time()-t:.1f}s  {sorted(res2)[:8]}", flush=True)
    if not res2:
        print("  *** ALL ROOTS 0 mod p ***")
        pickle.dump({'valp': env.valp[:], 'forced': dict(env.forced)},
                    open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_modp_sol.pkl','wb'))
        print("  saved mod-p solution")
        break
    if incons:
        print("  inconsistent - stop"); break
