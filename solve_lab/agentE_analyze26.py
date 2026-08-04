#!/usr/bin/env python3
"""Analyze the 26 breaks after control-zeroing. Classify each (square? vars wire/free/msgbit).
Find free-input handles in cones. Key: which handles move a break WITHOUT disturbing loads/core."""
import json, re, sys
from collections import defaultdict, deque
from agentE_common import build_wire, load_gates, load_lines, p, NVARS, VAR, CORE
sys.setrecursionlimit(1000000)
wire, find2, A = build_wire()
wire_set = set(wire)
gates = load_gates()
lines = load_lines()
eqcode = [compile(VAR.sub(r'v[\1]', L.rsplit('=', 1)[0]), '<e>', 'eval') for L in lines]
eqvars = [set(int(m) for m in VAR.findall(L)) for L in lines]
gate_out = set(t for t, _, _ in gates); defn = {}
for t, rhs, vids in gates:
    if t not in defn: defn[t] = vids
freeinp = set(v for v in range(NVARS) if v not in gate_out)
msgbit = set(t for t, rhs, vids in gates if rhs.strip() in ('0', '1'))

CORE_BRK = [11854, 29437, 32916]
NC_BRK = [3408,3841,4134,4526,5069,7276,15440,15724,15927,21600,22139,22825,27289,27999,28718,29305,31134,31269,32463,33195,36387,36390,38888]
ALL = CORE_BRK + NC_BRK

anc_cache = {}
def anc(v):
    if v in anc_cache: return anc_cache[v]
    if v in freeinp or v not in defn: r = {v} if v in freeinp else set(); anc_cache[v] = r; return r
    anc_cache[v] = set(); s = set()
    for u in defn[v]: s |= anc(u)
    anc_cache[v] = s; return s

# raw structure: is it a literal square A*A?
def is_square(i):
    L = lines[i].rsplit('=', 1)[0]
    # crude: strip leading const*, check top-level '(X)*(X)'
    return '*' in L

print("=== 26 break classification ===")
cone_union = set()
for i in ALL:
    vs = eqvars[i]
    w = len(vs & wire_set); fr = sorted(vs & freeinp); mb = len(vs & msgbit)
    cone = set()
    for v in vs: cone |= anc(v)
    cone -= wire_set
    cone_union |= cone
    tag = 'CORE' if i in CORE_BRK else ''
    print(f" eq {i}{tag}: {len(vs)}v wire{w} free{len(fr)}:{fr[:4]} msg{mb} | cone_free={len(cone)}")

print(f"\nunion of free-input cones (candidate handles): {len(cone_union)}")
# which of these handles are in the LOAD cones (would disturb loads if changed)?
load_cone = set()
for lv in (11150, 25739, 37758, 35389, 6671):
    load_cone |= anc(lv)
clean = cone_union - load_cone - {14853, 16742, 12186}
print(f"handles NOT in load cone (safe, won't disturb loads): {len(clean)}")
print(f"handles IN load cone (would disturb loads/core): {len(cone_union & load_cone)}")
# show a sample of clean handles with footprint
occ = defaultdict(int)
for L in lines:
    for m in set(int(x) for x in VAR.findall(L)): occ[m] += 1
sample = sorted(clean, key=lambda h: occ[h])[:25]
print("clean handles (sorted by footprint): " + ", ".join(f"x_{h}({occ[h]})" for h in sample))
