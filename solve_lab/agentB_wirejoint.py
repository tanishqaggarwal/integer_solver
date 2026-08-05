#!/usr/bin/env python3
"""From wire1correct base: fix ONLY the 3 true core products (wire 5101,32017,26789 + handles
30317,2936,5146). Joint linear solve over (all other wire + all other handles) to heal the 13
unpacking + identity-ripple, keeping everything else satisfied (exclude 15 M2-core which need the
L2 fix). Check consistency, apply, and re-evaluate the ACTUAL failing (tests bilinear 2nd-order)."""
import json, pickle, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire)
sol = {int(k[2:]): int(v) for k, v in json.load(open('wire1correct.json')).items()}
FIXW = {5101, 32017, 26789}      # wire members in the 3 core products
FIXH = {30317, 2936, 5146}       # core quotient handles
M2CORE = set([2071,7123,7469,13660,15299,16622,17726,21382,22093,25480,28653,31061,32894,34517,34892])
UNPACK = set([8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666])

# base: load wire1correct fully; wire forced=1, handles from json
env.forced = {v: (sol.get(v, 1) % p) for v in wire}   # wire members at their json (=1) values
env.jac_free = wireset - FIXW
env.set_free({v: sol.get(v, 0) for v in env.freeset})
# core handles already in freevals from json; verify base
res0 = env.all_root_residuals()
print(f"[wj] base failing={len(res0)}  (unpack={len(set(res0)&UNPACK)}, M2core={len(set(res0)&M2CORE)}, other={len(set(res0)-UNPACK-M2CORE)})")

allowed = (wireset - FIXW) | (set(env.freeinp) - FIXH)
env.tangent_linear()
rows = []
targets = 0
for i in range(len(env.root_poly)):
    if i in M2CORE:   # leave M2-core out (handled via L2 later)
        continue
    g = env.root_grad(i)
    if not g: continue
    gr = {c: v for c, v in g.items() if c in allowed}
    r = env.root_val(i)
    if not gr:
        if r % p: pass  # unfixable with this block
        continue
    rows.append((i, gr, (-r) % p))
coldeg = defaultdict(int)
for _, rd, _ in rows:
    for c in rd: coldeg[c] += 1
pivots = {}; piv_order = []; incons = 0; ilist=[]
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
    if not rd:
        if rhs % p: incons += 1; ilist.append(i)
        continue
    pc = min(rd, key=lambda c: coldeg.get(c,0)); inv = pow(rd[pc],p-2,p)
    pivots[pc] = ({c:(v*inv)%p for c,v in rd.items()}, (rhs*inv)%p); piv_order.append(pc)
print(f"[wj] joint linear system: rank={len(pivots)}, INCONS={incons} {ilist[:8]}")
delta = {}
for pc in reversed(piv_order):
    prow, prhs = pivots[pc]; s = prhs
    for c, v in prow.items():
        if c != pc:
            dv = delta.get(c,0)
            if dv: s = (s - v*dv) % p
    if s: delta[pc] = s
# apply
for c, d in delta.items():
    if c in wireset: env.forced[c] = (env.forced[c] + d) % p
    else: env.valp[c] = (env.valp[c] + d) % p
env.forward()
res2 = env.all_root_residuals()
print(f"[wj] after applying delta (|delta|={len(delta)}): ACTUAL failing={len(res2)}")
print(f"     unpack={len(set(res2)&UNPACK)}, M2core={len(set(res2)&M2CORE)}, other={sorted(set(res2)-UNPACK-M2CORE)[:20]}")
