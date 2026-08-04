#!/usr/bin/env python3
"""Two-phase from agentA+wire=1 (core handles fixed at 0):
 Phase 1: one wire+handle step heals the 13 unpacking (wire moves).
 Phase 2: FIX the wire, Newton on handles only (linear, no wire*handle 2nd-order) to heal the
          remaining 12 ripple failures.  If all -> 0 mod p, save for Dixon lift."""
import json, pickle, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire)
COREH = {30317, 2936, 5146}
sol = {int(k[2:]): int(v) for k, v in json.load(open('best_agentA_39021.json')).items()}
freevals = {v: sol.get(v, 0) for v in env.freeset}

env.forced = {v: (s % p) for v, s in wire.items()}
for h in COREH: env.forced[h] = 0

def step(jfree):
    env.jac_free = jfree
    env.tangent_linear()
    res = env.all_root_residuals()
    rows = []
    for i in range(len(env.root_poly)):
        g = env.root_grad(i)
        if g: rows.append((i, g))
    coldeg = defaultdict(int)
    for _, rd in rows:
        for c in rd: coldeg[c] += 1
    pivots = {}; piv_order = []; incons = 0; ilist = []
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
    for pc in reversed(piv_order):
        prow, prhs = pivots[pc]; s = prhs
        for c, v in prow.items():
            if c != pc:
                dv = delta.get(c,0)
                if dv: s = (s - v*dv) % p
        if s: delta[pc] = s
    return len(res), incons, ilist, delta

def apply(delta):
    for c, d in delta.items():
        if c in wireset: env.forced[c] = (env.forced[c] + d) % p
        else: env.valp[c] = (env.valp[c] + d) % p
    env.forward()

env.set_free(freevals)
print(f"[A3] start fail={len(env.all_root_residuals())}")
# Phase 1: wire+handles, one step (heals unpacking)
nfail, incons, ilist, delta = step(wireset)
apply(delta)
print(f"[A3] phase1 (wire+handles): {nfail}->{len(env.all_root_residuals())} incons={incons} |delta|={len(delta)}")

# Phase 2: FIX wire (jac_free empty -> only free-input handles move), Newton on handles
for it in range(30):
    t = time.time()
    nfail, incons, ilist, delta = step(set())   # no wire columns
    apply(delta)
    res2 = env.all_root_residuals()
    print(f"  H-it {it}: {nfail}->{len(res2)} incons={incons} {ilist[:5]} |delta|={len(delta)} {time.time()-t:.1f}s", flush=True)
    if not res2:
        print("  *** ALL 39033 ROOTS 0 mod p ***")
        pickle.dump({'valp': env.valp[:], 'forced': dict(env.forced), 'freevals': freevals,
                     'wireset': list(wireset)},
                    open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_modp_solA.pkl','wb'))
        print("  saved mod-p solution"); break
    if incons and it >= 2:
        print("  handle-only inconsistent - stop"); break
