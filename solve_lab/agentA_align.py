#!/usr/bin/env python3
"""Alignment solve: find live additive free knobs for x_19083 and x_23927, set them so
x_19083=x_24908 and x_23927=x_14853, then x_16742=x_24908, x_12186=x_14853. This keeps
ALL equality checks (x_16742=x_19083, x_12186=x_23927) AND zeroes x_3558,x_29322."""
import json, sys
from agentA_harness import (p, load_solution, forward, gates, definer, order, freeinp,
                            backward_cone, eqcode, eqvars, NEQ, NVARS)
base = load_solution('best/new_instance_partial_39013.json'); forward(base)
def count(v):
    ns = {'__builtins__': {}, 'v': v}
    return set(i for i in range(NEQ) if eval(eqcode[i], ns) != 0)
F0 = count(base)
gvids = {t: gates[definer[t]][2] for t in order}

def free_cone(root):
    seen=set(); st=[root]
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        for w in gvids.get(u, ()):
            if w not in seen: st.append(w)
    return [u for u in seen if u in freeinp]

def find_additive_knob(gate, avoid):
    """find a free input w with d(gate)/dw == +-1 (live additive), not in avoid."""
    v = base[:]; forward(v); g0 = v[gate]
    for w in free_cone(gate):
        if w in avoid: continue
        old = v[w]; v[w] = old + 1; forward(v); d = v[gate] - g0; v[w] = old
        if d in (1, -1):
            forward(v)
            return w, d
    forward(base)
    return None, None

# S,T critical + quotient protect
crit = {16742, 14853, 12186, 24908}
for r in [35389, 6671]:
    _, fr = backward_cone(r); crit |= set(fr)

k19, d19 = find_additive_knob(19083, crit)
k23, d23 = find_additive_knob(23927, crit)
print(f"x_19083 additive knob: x_{k19} (deriv {d19})")
print(f"x_23927 additive knob: x_{k23} (deriv {d23})")

# Build aligned state
def build():
    v = base[:]
    forward(v)  # baseline gates
    # target: x_19083 = x_24908 ; x_23927 = x_14853
    if k19 is not None:
        v[k19] = base[k19] + d19 * (v[24908] - v[19083])
    if k23 is not None:
        v[k23] = base[k23] + d23 * (v[14853] - v[23927])
    forward(v)
    # now set the free pins
    v[16742] = v[24908]     # x_3558 = 0, and x_16742 = x_19083(=x_24908)
    v[12186] = v[14853]     # x_29322 = 0, and x_12186 = x_23927(=x_14853)
    forward(v)
    if v[11150] % p == 0: v[30317] = -(v[11150])//p
    if (537773*v[37758]) % p == 0: v[2936] = (537773*v[37758])//p
    if v[25739] % (6672769*p) == 0: v[5146] = v[25739]//(6672769*p)
    return v

v = build()
print(f"\nafter align+pin:")
print(f"  x_19083==x_24908? {v[19083]==v[24908]}   x_16742==x_19083? {v[16742]==v[19083]}")
print(f"  x_23927==x_14853? {v[23927]==v[14853]}   x_12186==x_23927? {v[12186]==v[23927]}")
print(f"  x_3558 mod p={v[3558]%p}  x_29322 mod p={v[29322]%p}  S,T0={v[35389]%p==0},{v[6671]%p==0}")
F = count(v)
broke = sorted(F - F0); fixed = sorted(F0 - F)
print(f"  satisfied {NEQ-len(F)}/{NEQ} ({len(F)} fail)")
print(f"  broke {len(broke)}: {broke}")
print(f"  still-failing: {sorted(F)}")
if NEQ - len(F) > 39016:
    json.dump({f"x_{i}": v[i] for i in range(NVARS) if v[i] != 0}, open(f'best_agentA_{NEQ-len(F)}.json', 'w'))
    print(f"  SAVED best_agentA_{NEQ-len(F)}.json")
