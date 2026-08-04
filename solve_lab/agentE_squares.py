#!/usr/bin/env python3
"""Analyze the 17 decouple-break equations. For each, classify its DIRECT variables (wire / free /
message-bit / gate-output) and, via forward-eval perturbation, find free-input handles in the cone
that can re-zero it. Uses the (Q,g) trick for squares: eq=Q^2, perturb handle to solve Q->0."""
import json, re, sys
from collections import defaultdict, deque
from math import isqrt
from agentE_common import build_wire, load_gates, load_lines, p, NVARS, VAR, CORE
sys.setrecursionlimit(1000000)
M2MOD = 6672769
wire, find2, A = build_wire()
wire_set = set(wire)
gates = load_gates()
lines = load_lines()
eqcode = [compile(VAR.sub(r'v[\1]', L.rsplit('=', 1)[0]), '<e>', 'eval') for L in lines]
eqvars = [set(int(m) for m in VAR.findall(L)) for L in lines]
gate_out = set(t for t, _, _ in gates)
freeinp = set(v for v in range(NVARS) if v not in gate_out)
defn = {}
for t, rhs, vids in gates:
    if t not in defn: defn[t] = vids
best = {int(k[2:]): v for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}

BREAKS = [765,1658,10743,16616,19232,23663,25003,25455,26999,28742,30834,32082,36331,36442,36950,37697,38953]

# free-input ancestors (memoized)
anc_cache = {}
def anc(v):
    if v in anc_cache: return anc_cache[v]
    if v in freeinp or v not in defn: r = {v} if v in freeinp else set(); anc_cache[v] = r; return r
    anc_cache[v] = set(); s = set()
    for u in defn[v]: s |= anc(u)
    anc_cache[v] = s; return s

# classify message bits: gate outputs with constant 0/1 rhs
msgbit = set()
for t, rhs, vids in gates:
    if rhs.strip() in ('0', '1'): msgbit.add(t)

print("=== break-equation direct-variable classification ===")
for i in BREAKS:
    vs = eqvars[i]
    w = sorted(vs & wire_set); fr = sorted(vs & freeinp); mb = sorted(vs & msgbit)
    other = sorted(vs - wire_set - freeinp - msgbit)
    core3 = sorted(vs & {5101, 32017, 26789})
    print(f" eq {i}: {len(vs)} vars | wire {len(w)} | free {len(fr)}:{fr[:6]} | msgbit {len(mb)} | core3 {core3} | othergate {len(other)}")

# For each break, union free-input cone, find handles that move the eq. Use (Q,g) square trick.
# Build a minimal forward-eval limited to needed gates. Simpler: full forward-eval harness (reuse pattern).
fixed = {5101, 32017, 26789, 30317, 5146, 2936}
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
for m in (5101, 32017, 26789): val[m] = 1
ns = {'__builtins__': {}, 'v': val}
def forward():
    ns['v'] = val
    for k, t in enumerate(order): val[t] = eval(gcode[k], ns)
forward()
def evi(i): ns['v'] = val; return eval(eqcode[i], ns)

print("\n=== handle probe per break (free-input cone) ===")
for i in BREAKS:
    cone = set()
    for v in eqvars[i]: cone |= anc(v)
    cone -= wire_set
    R = evi(i)
    # Q trick to get |Q|,g for a handle
    found = []
    for h in sorted(cone):
        v0 = val[h]
        e0 = R
        val[h] = v0+1; forward(); e1 = evi(i)
        val[h] = v0+2; forward(); e2 = evi(i)
        val[h] = v0; forward()
        d1 = e1-e0; d2 = e2-e1
        gg2 = (d2-d1)  # = 2 g^2 for a pure square; but eq may be k*Q^2 -> 2 k g^2
        if d1 == 0 and d2 == 0: continue
        found.append((h, d1, d2))
        if len(found) >= 6: break
    print(f" eq {i}: R bits={R.bit_length()}, cone free={len(cone)}, movers(sample)={[(h,) for h,_,_ in found][:6]}")
