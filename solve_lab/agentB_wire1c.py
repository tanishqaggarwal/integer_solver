#!/usr/bin/env python3
"""Investigate the 13 unpacking equations at wire=1: mod-p residuals, trace the gate terms to
their free-input ancestors, and test whether setting those ancestors can heal them."""
import json, pickle, random
from collections import defaultdict, deque
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
env.forced = {v: (s % p) for v, s in wire.items()}
env.set_from_solution(best)
UNPACK = [8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666]
print("[c] mod-p residuals of the 13 unpacking eqs at wire=1:")
for i in UNPACK:
    r = env.root_val(i)
    # signed small representation
    rs = r if r < p//2 else r - p
    print(f"   eq {i}: resid mod p = {rs}")

# build gate input map (var -> its gate's input vars) from gate_poly
gate_inputs = {}
for t, pol in env.gate_poly:
    vs = set()
    for m in pol:
        vs |= set(m)
    gate_inputs[t] = vs
freeset = env.freeset

def free_ancestors(v, maxdepth=40):
    seen = set(); anc = set(); stack = [(v, 0)]
    while stack:
        x, d = stack.pop()
        if x in seen: continue
        seen.add(x)
        if x in freeset: anc.add(x); continue
        if x in wire: continue  # forced constant
        for u in gate_inputs.get(x, ()):
            stack.append((u, d+1))
    return anc

# for each unpacking eq, collect the gate vars in its root and their free ancestors
for i in UNPACK:
    rp = env.root_poly[i]
    gatevars = set()
    for m in rp:
        for v in m:
            if v not in wire and v not in freeset:
                gatevars.add(v)
    fanc = set()
    for gv in gatevars: fanc |= free_ancestors(gv)
    print(f"[c] eq {i}: {len(gatevars)} gate-vars, {len(fanc)} free ancestors {sorted(fanc)[:12]}")

# Focus eq 8429: try toggling each free ancestor and see effect on residual
i0 = 8429
rp = env.root_poly[i0]
gatevars = set(v for m in rp for v in m if v not in wire and v not in freeset)
fanc = set()
for gv in gatevars: fanc |= free_ancestors(gv)
print(f"\n[c] eq {i0}: free ancestors = {sorted(fanc)}")
base_val = {v: env.valp[v] for v in freeset}
r0 = env.root_val(i0)
moved = []
for h in sorted(fanc):
    old = env.valp[h]
    for testv in [(old+1)%p, (old+7)%p, 1, 0, (p-1)]:
        env.valp[h] = testv; env.forward()
        if env.root_val(i0) != r0:
            moved.append(h); break
    env.valp[h] = old; env.forward()
print(f"[c] eq {i0}: free ancestors that CHANGE the residual: {moved}")

# what are the gate vars and their current values / definitions
print(f"[c] eq {i0} gate vars and values:")
for gv in sorted(gatevars):
    ins = gate_inputs.get(gv, set())
    inkind = {('W' if v in wire else 'F' if v in freeset else 'G'): v for v in ins}
    print(f"   x_{gv}={env.valp[gv]%p if env.valp[gv]%p<10**6 else '...'}  inputs={sorted(ins)} ({inkind})")
