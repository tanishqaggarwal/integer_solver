#!/usr/bin/env python3
"""Solve the 228 pure-wire equations for the 220 wire members (Newton over wire only). Decisive:
is there a wire assignment satisfying ALL pure-wire eqs with the 3 core-product members != 0?
Try (A) core pinned=1, solve rest; (B) all wire free from wire=1."""
import json, pickle
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire); freeset = env.freeset
COREW = {5101, 32017, 26789}
purewire = [i for i, vs in enumerate(env.eqvars) if (vs & wireset) and not (vs & freeset)]
print(f"[pw] pure-wire equations: {len(purewire)}")

env.forced = {v: (s % p) for v, s in wire.items()}
env.set_free({v: best.get(v, 0) for v in env.freeset})

def solve_wire(jfree, tag):
    env.jac_free = jfree
    for it in range(30):
        env.tangent_linear()
        rows = []; nfail = 0
        for i in purewire:
            r = env.root_val(i)
            if r: nfail += 1
            g = env.root_grad(i)
            gr = {c: v for c, v in g.items() if c in jfree}
            if not gr:
                continue
            rows.append((gr, (-r) % p))
        if nfail == 0:
            print(f"[{tag}] it {it}: pure-wire SOLVED (0 fail). core wire: "
                  f"x5101={env.valp[5101]} x32017={env.valp[32017]} x26789={env.valp[26789]}")
            return True
        coldeg = defaultdict(int)
        for rd, _ in rows:
            for c in rd: coldeg[c] += 1
        pivots = {}; piv_order = []; incons = 0
        for k in sorted(range(len(rows)), key=lambda k: len(rows[k][0])):
            rd = dict(rows[k][0]); rhs = rows[k][1]
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
        for c, d in delta.items(): env.forced[c] = (env.forced[c] + d) % p
        env.forward()
        if it < 3 or it % 5 == 0:
            print(f"[{tag}] it {it}: pure-wire fail={nfail} incons={incons} |d|={len(delta)}", flush=True)
        if incons and it >= 2:
            print(f"[{tag}] it {it}: INCONSISTENT (fail={nfail}, incons={incons})")
            return False
    return None

# (A) core pinned=1
env.forced = {v: (s % p) for v, s in wire.items()}
env.set_free({v: best.get(v, 0) for v in env.freeset})
print("\n=== (A) core wire pinned=1, solve other 217 for pure-wire ===")
solve_wire(wireset - COREW, "A")
print(f"[A] final core wire: x5101={env.valp[5101]} x32017={env.valp[32017]} x26789={env.valp[26789]}")

# check: does the pure-wire homogeneous system have nullspace allowing core != 0?
# rank of pure-wire Jacobian (linear part at wire=1) over all 220
env.forced = {v: (s % p) for v, s in wire.items()}
env.set_free({v: best.get(v, 0) for v in env.freeset})
env.jac_free = wireset
env.tangent_linear()
rows = []
for i in purewire:
    g = env.root_grad(i)
    gr = {c: v for c, v in g.items() if c in wireset}
    if gr: rows.append(gr)
# rank
coldeg = defaultdict(int)
for rd in rows:
    for c in rd: coldeg[c] += 1
pivots = set()
piv = {}
for rd0 in sorted(rows, key=len):
    rd = dict(rd0)
    for c in list(rd):
        if c in piv:
            f = rd[c]; prow = piv[c]
            for cc, v in prow.items():
                nv = (rd.get(cc,0)-f*v) % p
                if nv: rd[cc]=nv
                elif cc in rd: del rd[cc]
    rd = {c:v for c,v in rd.items() if v}
    if rd:
        pc = min(rd); inv=pow(rd[pc],p-2,p)
        piv[pc]={c:(v*inv)%p for c,v in rd.items()}
print(f"\n[pw] pure-wire Jacobian rank over 220 wire = {len(piv)}; nullity = {220-len(piv)}")
print(f"[pw] => wire solution space dimension = {220-len(piv)} (if 0, wire is rigidly forced)")
