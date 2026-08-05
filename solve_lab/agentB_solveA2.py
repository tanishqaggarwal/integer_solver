#!/usr/bin/env python3
"""From agentA(core solved, loads=0)+wire=1: fix the 3 core handles at 0 (kills core 2nd-order),
let all 220 wire + all other handles move. Check consistency, single-step heal, then iterate."""
import json, pickle, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire)
COREH = {30317, 2936, 5146}
sol = {int(k[2:]): int(v) for k, v in json.load(open('best_agentA_39021.json')).items()}
freevals = {v: sol.get(v, 0) for v in env.freeset}
UNPACK = set([8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666])
CORE = set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])

# forced: wire=1 (jac_free), core handles=0 (constant, not jac_free)
env.forced = {v: (s % p) for v, s in wire.items()}
for h in COREH: env.forced[h] = 0
env.jac_free = wireset          # core handles NOT in jac_free -> gradient 0 -> stay fixed
env.set_free(freevals)
res = env.all_root_residuals()
fc = sorted(set(res)&CORE); fu = sorted(set(res)&UNPACK); fo = sorted(set(res)-CORE-UNPACK)
print(f"[A2] wire=1, coreH=0: fail={len(res)} core={len(fc)} unpack={len(fu)} other={len(fo)} {fo}")

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

for it in range(30):
    t = time.time()
    nfail, incons, ilist, delta = step()
    for c, d in delta.items():
        if c in wireset: env.forced[c] = (env.forced[c] + d) % p
        else: env.valp[c] = (env.valp[c] + d) % p
    env.forward()
    res2 = env.all_root_residuals()
    print(f"  it {it}: fail {nfail}->{len(res2)}, incons={incons} {ilist[:4]}, |delta|={len(delta)}, {time.time()-t:.1f}s", flush=True)
    if not res2:
        print("  *** ALL 39033 ROOTS 0 mod p ***")
        pickle.dump({'valp': env.valp[:], 'forced': dict(env.forced), 'freevals': freevals, 'wireset': list(wireset)},
                    open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_modp_solA.pkl','wb'))
        print("  saved mod-p solution"); break
    if incons and it >= 2:
        print("  inconsistent - stop"); break
