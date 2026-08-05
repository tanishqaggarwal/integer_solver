#!/usr/bin/env python3
"""Wire SPLIT: solve the 13 unpacking (linear in wire) for a wire pattern with the 3 core-product
members kept =1, then FIX the wire and handle-solve the ripple + core (linear, wire fixed)."""
import json, pickle, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire)
CORE = set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
UNPACK = set([8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666])
COREW = {5101, 32017, 26789}; COREH = {30317, 2936, 5146}
handles = set(env.freeinp)
inv6672769 = pow(6672769, p-2, p)
env.forced = {v: (s % p) for v, s in wire.items()}
env.jac_free = wireset
env.set_free({v: best.get(v, 0) for v in env.freeset})

# --- Step 1: solve the 13 unpacking for wire values, core members fixed=1 ---
# unpacking roots are LINEAR in wire (given handles fixed). Solve A w = b over wire (excl COREW).
env.tangent_linear()
solve_w = wireset - COREW
rows = []
for i in UNPACK:
    g = env.root_grad(i)
    gr = {c: v for c, v in g.items() if c in solve_w}
    r = env.root_val(i)
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
print(f"[split] unpacking wire-solve: rank={len(pivots)}, incons={incons}, |delta|={len(delta)}")
for c, d in delta.items(): env.forced[c] = (env.forced[c] + d) % p
env.forward()
res = env.all_root_residuals()
print(f"[split] after wire pattern: fail={len(res)} unpk={len(set(res)&UNPACK)} core={len(set(res)&CORE)} "
      f"other={len(set(res)-set(UNPACK)-CORE)}; coreW: x5101={env.valp[5101]} x32017={env.valp[32017]} x26789={env.valp[26789]}")

# --- Step 2: FIX wire, handle-solve ripple + core ---
def set_core_handles():
    vp = env.valp
    vp[30317] = (-vp[11150] * pow(vp[5101], p-2, p)) % p
    vp[2936] = (537773 * vp[37758] % p) * pow(vp[26789], p-2, p) % p
    vp[5146] = vp[25739] * inv6672769 % p * pow(vp[32017], p-2, p) % p
    env.forward()
env.jac_free = set()  # wire fixed
for it in range(15):
    set_core_handles()
    env.tangent_linear()
    resd = env.all_root_residuals()
    rows = []
    for i in range(len(env.root_poly)):
        if i in CORE: continue
        g = env.root_grad(i)
        if not g: continue
        gr = {c: v for c, v in g.items() if c in handles}
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
    for c, d in delta.items(): env.valp[c] = (env.valp[c] + d) % p
    env.forward(); set_core_handles()
    res2 = env.all_root_residuals()
    print(f"  s2-it {it}: fail={len(res2)} core={len(set(res2)&CORE)} oth={len(set(res2)-set(UNPACK)-CORE)} incons={incons} |d|={len(delta)}", flush=True)
    if not res2:
        print("  *** ALL 39033 ROOTS 0 mod p ***")
        pickle.dump({'valp': env.valp[:], 'forced': dict(env.forced)},
                    open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_modp_sol.pkl','wb'))
        print("  SAVED"); break
    if incons and it >= 2: print("  plateau"); break
