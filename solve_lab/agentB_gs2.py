#!/usr/bin/env python3
"""Alternating block GS, wire=1, 39013 base. H-step: handles solve core+pair-breaks. W-step: only
MOVABLE wire (excl core16) heals unpacking so the core stays solved. Categorize residual fails."""
import json, pickle, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire)
CORE16 = {3915, 5101, 10603, 16153, 17952, 21023, 22399, 23917, 26789, 30022, 31521, 32017, 32762, 35409, 38100, 38571}
handles = set(env.freeinp)
movable_wire = wireset - CORE16
UNPACK = set([8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666])
CORE = set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
env.forced = {v: (s % p) for v, s in wire.items()}
env.set_free({v: best.get(v, 0) for v in env.freeset})

def block_solve(allowed):
    env.jac_free = wireset
    env.tangent_linear()
    res = env.all_root_residuals()
    rows = []
    for i in range(len(env.root_poly)):
        g = env.root_grad(i)
        if not g: continue
        gr = {c: v for c, v in g.items() if c in allowed}
        rows.append((i, gr, (-res[i]) % p if i in res else 0))
    coldeg = defaultdict(int)
    for _, rd, _ in rows:
        for c in rd: coldeg[c] += 1
    pivots = {}; piv_order = []
    for k in sorted(range(len(rows)), key=lambda k: len(rows[k][1])):
        i, rd0, rhs = rows[k]; rd = dict(rd0)
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
        if not rd: continue
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
    return delta

def apply(delta):
    for c, d in delta.items():
        if c in wireset: env.forced[c] = (env.forced[c] + d) % p
        else: env.valp[c] = (env.valp[c] + d) % p
    env.forward()

def cats():
    res = env.all_root_residuals()
    return len(res), sorted(set(res)&CORE), sorted(set(res)&UNPACK), sorted(set(res)-CORE-UNPACK)

print(f"[gs2] start: fail={cats()[0]}")
for rnd in range(50):
    t = time.time()
    apply(block_solve(handles))
    apply(block_solve(movable_wire))
    n, fc, fu, fo = cats()
    print(f"  rnd {rnd}: fail={n} core={len(fc)} unpack={len(fu)} other={len(fo)} {fo[:12]}  {time.time()-t:.1f}s", flush=True)
    if n == 0:
        print("  *** ALL 39033 ROOTS 0 mod p ***")
        pickle.dump({'valp': env.valp[:], 'forced': dict(env.forced),
                     'freevals': {v: env.valp[v] for v in env.freeset}, 'wireset': list(wireset)},
                    open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_modp_sol.pkl','wb'))
        print("  saved mod-p solution"); break
    if rnd > 5 and n == prev:
        print("  plateau"); break
    prev = n
