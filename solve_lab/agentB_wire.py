#!/usr/bin/env python3
"""Compute the signed identity wire for the NEW instance (union-find over 2-term identity atoms).
Find the class containing 38100/26064/5101/32017/26789, report members+signs, and its current
forced value. Save wire dict for the wire=1 re-analysis."""
import json, pickle
from collections import Counter
from propagate import load_atoms, atom_vars
p = 2**256 - 2**32 - 977

A = load_atoms()
parent = {}; sgn = {}
def find(x):
    parent.setdefault(x, x); sgn.setdefault(x, 1)
    if parent[x] == x: return x, 1
    r, s = find(parent[x]); parent[x] = r; sgn[x] = sgn[x] * s
    return r, sgn[x]
def union(a, b, rel):
    ra, sa = find(a); rb, sb = find(b)
    if ra == rb: return
    parent[ra] = rb; sgn[ra] = rel * sb * sa
nid = 0
for poly in A:
    if len(poly) == 2:
        (m1, c1), (m2, c2) = list(poly.items())
        if len(m1) == 1 and len(m2) == 1 and abs(c1) == abs(c2):
            rel = -1 if (c1 > 0) == (c2 > 0) else 1
            union(m1[0], m2[0], rel); nid += 1
print(f"[wire] identity atoms: {nid}")
# class of 38100
r, _ = find(38100)
wire = {}
for x in list(parent):
    rr, s = find(x)
    if rr == r: wire[x] = s
print(f"[wire] class of 38100: size={len(wire)}; +:{sum(1 for s in wire.values() if s>0)} -:{sum(1 for s in wire.values() if s<0)}")
for probe in [38100, 26064, 5101, 32017, 26789]:
    print(f"   x_{probe} in wire: {probe in wire}, sign={wire.get(probe)}")

# single-var pins on wire members (what forces the wire value)
pins = {}
for poly in A:
    vs = atom_vars(poly)
    if len(vs) == 1:
        v = next(iter(vs)); c0 = poly.get((), 0); c1 = poly.get((v,), 0); c2 = poly.get((v, v), 0)
        if c2 == 0 and c1 != 0 and (-c0) % c1 == 0:
            pins[v] = (-c0) // c1
wire_pins = {v: pins[v] for v in wire if v in pins}
print(f"[wire] single-var pins on wire members: {len(wire_pins)}")
for v, val in list(wire_pins.items())[:6]:
    print(f"   x_{v} pinned to {val}  (=p? {val==p}, =p*sign? {val==p*wire[v]}, val/p={val//p if val%p==0 else 'n/a'})")

# how many gates vs free among wire members (need gate list)
gates_t = set()
with open('atoms/gates.jsonl') as f:
    for line in f:
        gates_t.add(json.loads(line)['t'])
wire_gate = [v for v in wire if v in gates_t]
wire_free = [v for v in wire if v not in gates_t]
print(f"[wire] wire members that are GATE outputs: {len(wire_gate)}; FREE inputs: {len(wire_free)}")
print(f"[wire] free wire members: {sorted(wire_free)[:20]}")

pickle.dump({'wire': wire, 'root': r}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl', 'wb'))
print("[wire] saved wire dict")
