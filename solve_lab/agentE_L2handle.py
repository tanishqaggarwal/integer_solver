#!/usr/bin/env python3
"""Find a handle to make 6672769 | L2 = x_25739. Trace the free-input cone of x_25739, L1, L3,
and the base gates S=x_35389, T=x_6671. Report each free input's coeff in L2 mod 6672769 and its
equation footprint. Also factor 6672769 and analyze L1,L2,L3 in terms of S,T."""
import json, re
from sympy import factorint
from collections import defaultdict, deque
from agentE_common import build_wire, load_best, load_lines, load_gates, p, NVARS, VAR

print(f"6672769 = {factorint(6672769)}")
print(f"537773  = {factorint(537773)}")

wire, find2, A = build_wire()
best = load_best()
gates = load_gates()
lines = load_lines()

# equation footprint
occ = defaultdict(int)
for L in lines:
    for m in set(int(x) for x in VAR.findall(L)): occ[m] += 1

# build gate DAG and topo order (like forward_construct)
gate_out = set(t for t, _, _ in gates)
freeinp = set(v for v in range(NVARS) if v not in gate_out)
definer = {}
for t, rhs, vids in gates:
    if t not in definer: definer[t] = (rhs, vids)

# backward cone: free-input ancestors of a var (memoized DFS)
import sys; sys.setrecursionlimit(1000000)
anc_cache = {}
def anc(v):
    if v in anc_cache: return anc_cache[v]
    if v in freeinp or v not in definer:
        anc_cache[v] = {v}; return {v}
    anc_cache[v] = set()  # guard against cycles
    s = set()
    for u in definer[v][1]: s |= anc(u)
    anc_cache[v] = s; return s

# S=x_35389, T=x_6671, loads
for name, v in [('S', 35389), ('T', 6671), ('L1', 11150), ('L2', 25739), ('L3', 37758)]:
    a = anc(v)
    print(f"\n{name}=x_{v}: {len(a)} free-input ancestors; best={best.get(v,0)} ({best.get(v,0).bit_length()} bits)")

# For L2's cone, compute each free input's linear coeff in x_25739 (numerical derivative), mod 6672769
L2cone = sorted(anc(25739))
print(f"\n=== L2 free-input cone: {len(L2cone)} inputs ===")
# evaluate x_25739 as function of val: need forward-eval of its subcircuit. Build topo order for cone.
# Simpler: numerical partial derivative by perturbing each free input by 1 and re-evaluating x_25739 via its gate subtree.
val = dict(best)
def eval_gate(v, override=None):
    # recursively eval gate value from current val + override
    if override and v in override: return override[v]
    if v in freeinp or v not in definer: return val.get(v, 0)
    rhs, vids = definer[v]
    expr = VAR.sub(lambda m: str(eval_gate(int(m.group(1)), override)), rhs)
    return eval(expr, {'__builtins__': {}})
# baseline
base25739 = best.get(25739, 0)
M = 6672769
handles = []
for fv in L2cone:
    ov = {fv: val.get(fv, 0) + 1}
    nv = eval_gate(25739, ov)
    coef = nv - base25739
    if coef != 0:
        handles.append((fv, coef % M, coef, occ[fv], fv in wire))
print(f"free inputs with nonzero d(L2)/d(input): {len(handles)}")
# sort by footprint (fewest equations = cleanest handle) and coprimality to M
handles.sort(key=lambda x: (x[3], abs(x[2])))
print("input -> (coef mod 6672769, coef_bits, #eqs, in_wire)  [sorted by footprint]")
for fv, cm, coef, oc, inw in handles[:40]:
    from math import gcd
    print(f"  x_{fv}: coefmod={cm}, coefbits={coef.bit_length()}, #eqs={oc}, wire={inw}, gcd(coef,M)={gcd(cm,M) if cm else M}")
