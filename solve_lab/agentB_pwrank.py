#!/usr/bin/env python3
"""Definitive pure-wire feasibility. Are the 228 pure-wire eqs linear in wire? Build the linear
system (at wire=1) over all 220 wire members, solve J*dw = -resid, report consistency and the
resulting core-member values. If healing unpacking forces core->0, the wire!=p (core!=0) branch is
mod-p infeasible."""
import json, pickle
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire)
COREW = {5101, 32017, 26789}
purewire = [i for i, vs in enumerate(env.eqvars) if (vs & wireset) and not (vs & env.freeset)]

# linearity check: any pure-wire eq with a wire*wire monomial?
nonlin = 0
for i in purewire:
    for m in env.root_poly[i]:
        wc = sum(1 for v in m if v in wireset)
        if wc >= 2: nonlin += 1; break
print(f"[pw] pure-wire eqs: {len(purewire)}; nonlinear-in-wire (wire*wire term): {nonlin}")

env.forced = {v: (s % p) for v, s in wire.items()}
env.jac_free = wireset
env.set_free({v: best.get(v, 0) for v in env.freeset})
env.tangent_linear()

# system J*dw = -resid over 220 wire cols
rows = []
for i in purewire:
    g = env.root_grad(i)
    gr = {c: v for c, v in g.items() if c in wireset}
    r = env.root_val(i)
    if gr or r % p:
        rows.append((i, gr, (-r) % p))

def rref(rows):
    coldeg = defaultdict(int)
    for _, rd, _ in rows:
        for c in rd: coldeg[c] += 1
    pivots = {}; piv_order = []; incons = 0; ilist = []
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
    return pivots, piv_order, incons, ilist

pivots, piv_order, incons, ilist = rref(rows)
rank = len(pivots)
print(f"[pw] pure-wire linear system: rank={rank}, nullity(over 220)={220-rank}, INCONS={incons} {ilist[:6]}")
# particular solution dw
dw = {}
for pc in reversed(piv_order):
    prow, prhs = pivots[pc]; s = prhs
    for c, v in prow.items():
        if c != pc:
            d = dw.get(c,0)
            if d: s = (s - v*d) % p
    if s: dw[pc] = s
# resulting core values = 1 + dw[core]
print(f"[pw] particular solution: core members after heal:")
for c in COREW:
    val = (1 + dw.get(c, 0)) % p
    print(f"    x_{c} = {val}  (zero mod p? {val % p == 0})")
# Is a core!=0 solution in the affine solution space? check if core cols are free (non-pivot)
freecols = [c for c in wireset if c not in pivots]
print(f"[pw] free (non-pivot) wire columns: {len(freecols)}; core members free? "
      f"{[c for c in COREW if c in freecols]}")
# If core members are PIVOTS determined to 0, infeasible. If free, we can set them !=0.
for c in COREW:
    if c in pivots:
        prow, prhs = pivots[c]
        # its value = prhs - sum(prow[other]*dw[other]) -- determined
        print(f"    x_{c} is a PIVOT (determined); forced value shown above")
    else:
        print(f"    x_{c} is FREE -> can set nonzero")
