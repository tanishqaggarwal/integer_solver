#!/usr/bin/env python3
"""PINNED-WIRE path: zero the control gates exactly. Set x_29322=0 (x_14853:=x_12186), x_3558=0
(x_16742:=x_24908) -> S=T=0 -> loads L1=L2=L3=0 -> core M1=M2=M3=0 with partners=0. No divisibility.
Cost: wiring ripple from moving x_14853,x_16742. Report exactly what breaks."""
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
eqvars = [set(int(m) for m in VAR.findall(L)) for L in lines]
gate_out = set(t for t, _, _ in gates)
freeinp = set(v for v in range(NVARS) if v not in gate_out)
best = {int(k[2:]): v for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}

PARTNERS = {30317, 5146, 2936}
CTRL = {14853, 16742}   # free inputs we set explicitly
fixed = PARTNERS | CTRL
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
qq = deque(gi for gi in range(len(gates)) if gu[gi] == 0)
while qq:
    gi = qq.popleft(); t, rhs, vids = gates[gi]
    if ready[t] or t in fixed: continue
    definer[t] = gi; order.append(t); ready[t] = True
    for gj in using[t]:
        gu[gj] -= 1
        if gu[gj] == 0: qq.append(gj)
gcode = [compile(VAR.sub(r'v[\1]', gates[definer[t]][1]), '<r>', 'eval') for t in order]
val = [0]*NVARS
for v in freeinp: val[v] = best.get(v, 0)
ns = {'__builtins__': {}, 'v': val}
def forward():
    ns['v'] = val
    for k, t in enumerate(order): val[t] = eval(gcode[k], ns)

# step 1: forward with base to read x_24908, x_12186
forward()
x12186 = val[12186]; x24908 = val[24908]
print(f"base: x_29322={val[29322]}, x_3558={val[3558]}, x_24908={x24908}, x_12186={x12186}")

# step 2: set x_14853 = x_12186 (x_29322=0), x_16742 = x_24908 (x_3558=0)
val[14853] = x12186
val[16742] = x24908
forward()
print(f"after ctrl-zero: x_29322={val[29322]}, x_3558={val[3558]}, S=x_35389={val[35389]}, T=x_6671={val[6671]}")
L1, L2, L3 = val[11150], val[25739], val[37758]
print(f"loads: L1={L1}, L2={L2}, L3={L3}")

# step 3: partners
val[30317] = -L1 // p if L1 % p == 0 else 0  # with L1=0 -> 0
val[2936] = (537773*L3) // p if (537773*L3) % p == 0 else 0
val[5146] = L2 // (M2MOD*p) if L2 % (M2MOD*p) == 0 else 0
# but with loads=0, all partners=0; M1=L1+p*x30317=0 etc. Actually set to make Mi=0:
# M1 = L1 + x_5101*x_30317; x_5101=p. -> x_30317 = -L1//p (L1=0 -> 0)
forward()

F = [i for i in range(len(lines)) if eval(eqcode[i], ns)]
c = [i for i in F if i in CORE]; nc = [i for i in F if i not in CORE]
print(f"\nRESULT: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail); core {len(c)}: {sorted(c)}")
print(f"noncore {len(nc)}: {sorted(nc)}")
# save for further repair
json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('agentE_stfix.json', 'w'))
print("saved agentE_stfix.json")
