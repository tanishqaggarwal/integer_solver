#!/usr/bin/env python3
"""DECOUPLE route: pin only x_5101,x_32017,x_26789 to a value (default 1), keep the rest of the wire
at p (forward-eval propagates x_26064=p through identity gates). See which equations break.
Expectation: only the alias squares (atoms 8035/8036/8037) + core (fixed later by partners)."""
import json, re, sys
from collections import defaultdict, deque
from agentE_common import build_wire, load_gates, load_lines, p, NVARS, VAR, CORE
sys.setrecursionlimit(1000000)
M2MOD = 6672769

wire, find2, A = build_wire()
wire_set = set(wire)
gates = load_gates()
lines = load_lines()
eqcode = [compile(VAR.sub(r'v[\1]', L.rsplit('=', 1)[0]), '<e>', 'eval') for L in lines]
gate_out = set(t for t, _, _ in gates)
freeinp = set(v for v in range(NVARS) if v not in gate_out)
best = {int(k[2:]): v for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}

CORE3 = {5101, 32017, 26789}
PARTNERS = {30317, 5146, 2936}
V3 = int(sys.argv[1]) if len(sys.argv) > 1 else 1

# build order excluding fixed = CORE3 (pinned) + PARTNERS (free, set explicitly). Wire NOT pinned (computed).
fixed = CORE3 | PARTNERS
ready = [False]*NVARS
for v in range(NVARS):
    if v not in gate_out or v in freeinp or v in fixed: ready[v] = True
gu = [0]*len(gates); using = defaultdict(list)
for gi, (t, rhs, vids) in enumerate(gates):
    u = 0
    for v in vids:
        if not ready[v]: u += 1
        using[v].append(gi)
    gu[gi] = u
definer = {}; order = []
q = deque(gi for gi in range(len(gates)) if gu[gi] == 0)
while q:
    gi = q.popleft(); t, rhs, vids = gates[gi]
    if ready[t] or t in fixed: continue
    definer[t] = gi; order.append(t); ready[t] = True
    for gj in using[t]:
        gu[gj] -= 1
        if gu[gj] == 0: q.append(gj)
gcode = [compile(VAR.sub(r'v[\1]', gates[definer[t]][1]), '<r>', 'eval') for t in order]

val = [0]*NVARS
for v in freeinp: val[v] = best.get(v, 0)
# pin the 3 core members
for m in CORE3: val[m] = V3
ns = {'__builtins__': {}, 'v': val}
def forward():
    ns['v'] = val
    for k, t in enumerate(order): val[t] = eval(gcode[k], ns)
forward()

# check wire propagated to p
wp = [w for w in wire if w not in CORE3 and abs(val[w]) != p]
print(f"wire members (non-core3) not at +-p: {len(wp)} (should be 0): {wp[:10]}")
print(f"x_26064={val[26064]==p and 'p' or val[26064]}, x_4376={val[4376]==p and 'p' or 'other'}, x_5101={val[5101]}")

L1, L2, L3 = val[11150], val[25739], val[37758]
print(f"L2 % 6672769 = {L2 % M2MOD}")

# fails WITHOUT partners set (core will fail); classify noncore
F = [i for i in range(len(lines)) if eval(eqcode[i], ns) != 0]
nc = [i for i in F if i not in CORE]; c = [i for i in F if i in CORE]
print(f"decouple V3={V3} (partners unset): {len(lines)-len(F)}/{len(lines)}, core {len(c)}, noncore {len(nc)}: {sorted(nc)}")

# the alias-square equations we expected
ALIAS = [765,1658,10743,16616,19232,23663,25455,28742,30834,32082,36331,36950]
print(f"expected alias eqs failing: {[i for i in ALIAS if i in F]}")
print(f"noncore NOT in alias list: {[i for i in nc if i not in ALIAS]}")
