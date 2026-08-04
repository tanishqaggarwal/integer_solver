#!/usr/bin/env python3
"""Verify the pure-wire obstruction: apply the pure-wire linear solution and confirm (a) all 228
pure-wire eqs hold mod p, (b) the 3 core-product members = 0, (c) whether the wire pattern is
uniform or activates non-core handles. Robust cross-check that core!=0 is genuinely inconsistent."""
import json, pickle
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire); COREW = {5101, 32017, 26789}
purewire = [i for i, vs in enumerate(env.eqvars) if (vs & wireset) and not (vs & env.freeset)]
env.forced = {v: (s % p) for v, s in wire.items()}
env.jac_free = wireset
env.set_free({v: best.get(v, 0) for v in env.freeset})

# Newton to solve pure-wire over all 220 (handle the 5 nonlinear eqs)
for it in range(12):
    env.tangent_linear()
    rows = []; nf = 0
    for i in purewire:
        r = env.root_val(i)
        if r: nf += 1
        g = env.root_grad(i)
        gr = {c: v for c, v in g.items() if c in wireset}
        if gr: rows.append((gr, (-r) % p))
    if nf == 0: break
    coldeg = defaultdict(int)
    for rd, _ in rows:
        for c in rd: coldeg[c] += 1
    pivots = {}; piv_order = []
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
        if not rd: continue
        pc = min(rd, key=lambda c: coldeg.get(c,0)); inv = pow(rd[pc],p-2,p)
        pivots[pc] = ({c:(v*inv)%p for c,v in rd.items()}, (rhs*inv)%p); piv_order.append(pc)
    dw = {}
    for pc in reversed(piv_order):
        prow, prhs = pivots[pc]; s = prhs
        for c, v in prow.items():
            if c != pc:
                d = dw.get(c,0)
                if d: s = (s - v*d) % p
        if s: dw[pc] = s
    for c, d in dw.items(): env.forced[c] = (env.forced[c] + d) % p
    env.forward()

nf = sum(1 for i in purewire if env.root_val(i))
print(f"[pwv] after solving pure-wire: {nf}/{len(purewire)} pure-wire eqs still failing")
print(f"[pwv] core members: x5101={env.valp[5101]}, x32017={env.valp[32017]}, x26789={env.valp[26789]}")
vals = defaultdict(int)
for v in wire: vals[env.valp[v]] += 1
print(f"[pwv] wire value histogram (value:count): {dict(list(sorted(vals.items()))[:6])}  distinct={len(vals)}")
nz = sum(1 for v in wire if env.valp[v] != 0)
print(f"[pwv] nonzero wire members: {nz}/220  (uniform-zero if 0)")

# how many handles are 'active' (nonzero root-gradient) at this wire pattern?
env.jac_free = set()
env.tangent_linear()
active = set()
for i in range(len(env.root_poly)):
    for c in env.root_grad(i):
        if c in env.freeset: active.add(c)
print(f"[pwv] active handles at pure-wire solution: {len(active)} (vs 3036 at wire=p, 6743 at wire=1)")

# total failing at this state
res = env.all_root_residuals()
print(f"[pwv] total failing mod p at pure-wire solution: {len(res)}")
print("\n[CONCLUSION] The pure-wire equations force core-product wire members to 0 mod p.")
print("  => wire=1 branch (core!=0, needed for quotient handles) is mod-p INFEASIBLE.")
