#!/usr/bin/env python3
"""Two-stage from 39013+wire=1.
Stage 1: linear solve over (wire + handles), SKIP core rows -> find wire pattern healing the 13
         unpacking + identity ripple (core wire members free).
Stage 2: FIX the wire; linear solve over handles only, targeting ALL 39033 (now linear -> exact).
Goal: all 39033 = 0 mod p, then save for Dixon lift."""
import json, pickle, time
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire)
CORE = set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
UNPACK = set([8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666])
handles = set(env.freeinp)
env.forced = {v: (s % p) for v, s in wire.items()}
env.set_free({v: best.get(v, 0) for v in env.freeset})

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

def apply(delta):
    for c, d in delta.items():
        if c in wireset: env.forced[c] = (env.forced[c] + d) % p
        else: env.valp[c] = (env.valp[c] + d) % p
    env.forward()

def cats():
    res = env.all_root_residuals()
    return len(res), len(set(res)&UNPACK), len(set(res)&CORE), len(set(res)-UNPACK-CORE)

# Stage 1
env.jac_free = wireset
print("[stage1] heal unpacking + ripple (skip core)")
for it in range(25):
    incons, delta = lin_solve(wireset | handles, CORE)
    apply(delta)
    n, u, c, o = cats()
    print(f"  s1-it {it}: fail={n} unpk={u} core={c} oth={o} incons={incons} |d|={len(delta)}", flush=True)
    if u == 0 and o == 0:
        print("  stage1 done: unpacking+ripple healed"); break
# check core wire nonzero
print(f"[stage1] core wire vals: x5101={env.valp[5101]}, x32017={env.valp[32017]}, x26789={env.valp[26789]}")

# Stage 2: fix wire, solve handles for everything (incl core)
env.jac_free = set()   # wire fixed
print("[stage2] fix wire, handle-only solve for all 39033")
for it in range(15):
    incons, delta = lin_solve(handles, set())
    apply(delta)
    n, u, c, o = cats()
    print(f"  s2-it {it}: fail={n} unpk={u} core={c} oth={o} incons={incons} |d|={len(delta)}", flush=True)
    if n == 0:
        print("  *** ALL 39033 ROOTS 0 mod p ***")
        pickle.dump({'valp': env.valp[:], 'forced': dict(env.forced)},
                    open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_modp_sol.pkl','wb'))
        print("  saved mod-p solution"); break
    if incons and it >= 2:
        print("  stage2 plateau"); break
