#!/usr/bin/env python3
"""Newton with closed-form core each step. From 39013+wire=1:
 (1) recompute core quotient handles in closed form so M1=M2=M3=0 mod p given current loads/wire;
 (2) linear solve over (all wire + non-core handles) to heal the 13 unpacking + identity ripple,
     excluding the 20 core rows (handled by step 1). Iterate. Goal: all 39033 = 0 mod p."""
import json, pickle, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire)
COREH = {30317, 2936, 5146}
CORE = set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
UNPACK = set([8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666])
inv6672769 = pow(6672769, p-2, p)
env.forced = {v: (s % p) for v, s in wire.items()}
env.jac_free = wireset
env.set_free({v: best.get(v, 0) for v in env.freeset})

def set_core_handles():
    vp = env.valp
    # x_30317 = -x_11150 * inv(x_5101);  x_2936 = 537773*x_37758*inv(x_26789)
    # x_5146 = x_25739 * inv(6672769*x_32017)
    vp[30317] = (-vp[11150] * pow(vp[5101], p-2, p)) % p
    vp[2936] = (537773 * vp[37758] % p) * pow(vp[26789], p-2, p) % p
    vp[5146] = vp[25739] * inv6672769 % p * pow(vp[32017], p-2, p) % p
    env.forward()

def lin_solve(allowed, skip):
    env.tangent_linear()
    res = env.all_root_residuals()
    rows = []
    for i in range(len(env.root_poly)):
        if i in skip: continue
        g = env.root_grad(i)
        if not g: continue
        gr = {c: v for c, v in g.items() if c in allowed}
        r = env.root_val(i)
        if not gr: continue
        rows.append((i, gr, (-r) % p))
    coldeg = defaultdict(int)
    for _, rd, _ in rows:
        for c in rd: coldeg[c] += 1
    pivots = {}; piv_order = []; incons = 0
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
    return incons, delta

allowed = (wireset | set(env.freeinp)) - COREH
for it in range(40):
    t = time.time()
    set_core_handles()
    nfail0 = len(env.all_root_residuals())
    incons, delta = lin_solve(allowed, CORE)
    for c, d in delta.items():
        if c in wireset: env.forced[c] = (env.forced[c] + d) % p
        else: env.valp[c] = (env.valp[c] + d) % p
    env.forward()
    set_core_handles()
    res2 = env.all_root_residuals()
    print(f"  it {it}: fail(after core-set)={nfail0} -> {len(res2)} incons={incons} |d|={len(delta)} "
          f"unpk={len(set(res2)&UNPACK)} core={len(set(res2)&CORE)} oth={len(set(res2)-UNPACK-CORE)} {time.time()-t:.1f}s", flush=True)
    if not res2:
        print("  *** ALL 39033 ROOTS 0 mod p ***")
        pickle.dump({'valp': env.valp[:], 'forced': dict(env.forced)},
                    open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_modp_sol.pkl','wb'))
        print("  saved"); break
