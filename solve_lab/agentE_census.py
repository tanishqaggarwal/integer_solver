#!/usr/bin/env python3
"""INTEGER handle census for the 13 codeword-unpacking equations. For each, list every variable,
its role (wire member / message-bit constant-gate / other gate / FREE INPUT), its value at best,
and its total equation footprint. A live integer quotient handle = a variable appearing in exactly
THIS equation (footprint 1) with a p-scale/invertible coefficient. Report whether any exists."""
import json, re
from collections import defaultdict
from agentE_common import build_wire, load_best, load_lines, load_gates, p, NVARS, VAR

wire, find2, A = build_wire()
wire_set = set(wire)
best = load_best()
lines = load_lines()
gates = load_gates()
gate_out = set(t for t, _, _ in gates)
freeinp = set(v for v in range(NVARS) if v not in gate_out)
# constant-gate (message bit / const) outputs
constgate = {}
for t, rhs, vids in gates:
    if not vids and rhs.strip().lstrip('-').isdigit():
        constgate[t] = int(rhs.strip())

occ = defaultdict(int)
for L in lines:
    for m in set(int(x) for x in VAR.findall(L)): occ[m] += 1

UNPACK13 = [8429, 11166, 11915, 12594, 23869, 25313, 26785, 31400, 32300, 36106, 36767, 37257, 37666]

def role(v):
    if v in wire_set: return f"WIRE(sgn{wire[v]:+d})"
    if v in constgate: return f"CONST={constgate[v]}"
    if v in freeinp: return "FREE"
    return "GATE"

print("=== INTEGER HANDLE CENSUS: 13 codeword-unpacking equations ===\n")
private_handles_total = 0
for i in UNPACK13:
    vs = sorted(set(int(m) for m in VAR.findall(lines[i].rsplit('=', 1)[0])))
    # footprint-1 vars = appear ONLY in this equation
    priv = [v for v in vs if occ[v] == 1]
    freev = [v for v in vs if v in freeinp]
    print(f"eq {i}: {len(vs)} variables")
    # tabulate roles
    rc = defaultdict(list)
    for v in vs: rc[role(v).split('(')[0].split('=')[0]].append(v)
    summary = ", ".join(f"{k}:{len(g)}" for k, g in sorted(rc.items()))
    print(f"   roles: {summary}")
    print(f"   FREE-input vars: {freev if freev else 'NONE'}")
    print(f"   footprint-1 (private) vars: {priv if priv else 'NONE'}")
    # min footprint among all vars
    mn = min((occ[v], v) for v in vs)
    print(f"   min footprint any var: x_{mn[1]} in {mn[0]} eqs (role {role(mn[1])})")
    private_handles_total += len(priv)
    print()

print("="*60)
print(f"TOTAL private (footprint-1) variables across all 13 unpackings: {private_handles_total}")
print(f"TOTAL free-input variables across all 13 unpackings: {sum(1 for i in UNPACK13 for v in set(int(m) for m in VAR.findall(lines[i])) if v in freeinp)}")
# The forcing atom x_26064 - p: footprint
print(f"\nx_26064 (the wire-pin forcing var): value==p? {best.get(26064)==p}, footprint {occ[26064]} eqs, role {role(26064)}")
print("Verdict: a live integer quotient handle would need footprint==1 with p-scale coeff.")
