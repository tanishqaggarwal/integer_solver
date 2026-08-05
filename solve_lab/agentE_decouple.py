#!/usr/bin/env python3
"""Decoupling analysis:
1. Per-equation V-sensitivity g_i(V) = eq_i with all wire=sign*V, else best. Fit poly in V.
2. Which wire members appear in each of the 13 breaking eqs + the 20 core.
3. Do the 3 core members x_5101,x_32017,x_26789 appear in the 13 unpackings?
4. For each core member, the identity (alias) atoms tying it to the wire and where they occur."""
import json, re
from collections import defaultdict
from agentE_common import build_wire, load_best, load_lines, p, NVARS, CORE, NONCORE13, VAR

wire, find2, A = build_wire()
best = load_best()
lines = load_lines()
wire_set = set(wire)
CORE3 = [5101, 32017, 26789]

# per-eq variable sets
eqvars = [set(int(m) for m in VAR.findall(L)) for L in lines]

# 1. V-sensitivity of the 13 + core: evaluate eq with wire=sign*V
def eval_eq_wireV(idx, Vval):
    val = defaultdict(int)
    for k, x in best.items(): val[k] = x
    for w, s in wire.items(): val[w] = s*Vval
    expr = VAR.sub(lambda m: str(val[int(m.group(1))]), lines[idx].rsplit('=', 1)[0])
    return eval(expr, {'__builtins__': {}})

print("=== V-sensitivity of the 13 breaking eqs (g_i(V) sampled) ===")
for idx in NONCORE13:
    vals = {Vv: eval_eq_wireV(idx, Vv) for Vv in [0, 1, 2, 3]}
    # fit: g(V) = c0 + c1 V + c2 V^2 + c3 V^3  from V=0,1,2,3
    g0, g1, g2, g3 = vals[0], vals[1], vals[2], vals[3]
    # finite differences
    d1 = g1-g0; d2 = g2-g1; d3 = g3-g2
    dd1 = d2-d1; dd2 = d3-d2
    ddd = dd2-dd1
    wm = sorted(eqvars[idx] & wire_set)
    has_core3 = [m for m in CORE3 if m in eqvars[idx]]
    print(f"  eq {idx}: g(0)bits={g0.bit_length()}, deg? c3={ddd//6 if ddd else 0} c2coef={dd1//2}, #wiremem={len(wm)}, has_core3={has_core3}")

print("\n=== Do the 3 core members appear in the 13 unpackings? ===")
for m in CORE3:
    inwhich = [idx for idx in NONCORE13 if m in eqvars[idx]]
    print(f"  x_{m}: in unpackings {inwhich}")

print("\n=== Wire members in each of the 13 unpackings ===")
for idx in NONCORE13:
    wm = sorted(eqvars[idx] & wire_set)
    print(f"  eq {idx}: {len(wm)} wire members: {wm[:20]}")

print("\n=== Identity (alias) atoms tying each core member to the wire ===")
# find 2-term identity atoms containing a core member
for m in CORE3:
    print(f"\n x_{m} (sign {wire[m]}):")
    for ai, poly in enumerate(A):
        vs = set()
        for mon in poly:
            for v in mon: vs.add(v)
        if m in vs and len(vs) == 2 and poly.get((), 0) == 0:
            nonconst = [mm for mm in poly if mm != ()]
            if all(len(mm) == 1 for mm in nonconst):
                other = (vs - {m}).pop()
                c1 = poly.get((m,), 0); c2 = poly.get((other,), 0)
                if abs(c1) == abs(c2):
                    inwire = other in wire_set
                    # how many equations use this atom
                    d = json.loads(open('atoms/poly_atoms.jsonl').readlines()[ai])
                    print(f"   atom {ai}: {poly.get((m,))}*x_{m} + {poly.get((other,))}*x_{other}  (x_{other} in wire={inwire}), used in {d.get('n_eq')} eqs: {d.get('eqs')}")
