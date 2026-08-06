#!/usr/bin/env python3
"""Core arithmetic: L1,L2,L3 values, divisibilities, the multiply gates for x_4007/x_29804/x_35605,
and the exact core equation structure. Determine what setting the partners requires."""
import json, re
from collections import defaultdict
from agentE_common import build_wire, load_best, load_lines, load_gates, p, NVARS, CORE, VAR

wire, find2, A = build_wire()
best = load_best()
gates = load_gates()
def g(v): return best.get(v, 0)

L1 = g(11150); L2 = g(25739); L3 = g(37758)
print("=== Loads ===")
print(f"L1 = x_11150 = {L1}  ({L1.bit_length()} bits)")
print(f"L2 = x_25739 = {L2}  ({L2.bit_length()} bits)")
print(f"L3 = x_37758 = {L3}  ({L3.bit_length()} bits)")
print(f"\nL1 % p = {L1 % p}   (p|L1? {L1%p==0})")
print(f"L3 % p = {L3 % p}   (p|L3? {L3%p==0})")
print(f"L2 % 6672769 = {L2 % 6672769}   (6672769|L2? {L2%6672769==0})")
print(f"L1 % 6672769 = {L1 % 6672769}")
print(f"L3 % 6672769 = {L3 % 6672769}")

print("\n=== Current core-related values ===")
for v in [15298, 4007, 29804, 35605, 30317, 5146, 2936, 5101, 32017, 26789]:
    print(f"  x_{v} = {g(v)}  ({g(v).bit_length()} bits)")

print("\n=== Multiply gates defining x_4007, x_29804, x_35605 ===")
for t, rhs, vids in gates:
    if t in (4007, 29804, 35605):
        print(f"  x_{t} = {rhs}   vids={vids}")
# also gates producing x_11150,x_25739,x_37758 (the loads)
print("\n=== Gates for loads x_11150,x_25739,x_37758 (may be sums) ===")
for t, rhs, vids in gates:
    if t in (11150, 25739, 37758):
        print(f"  x_{t} = {rhs}")

print("\n=== Core equation sanity: eq 22093 (=M2), 25539 (=7 M3), 2071 (=M1+5 M2) ===")
lines = load_lines()
def evline(idx, val):
    expr = VAR.sub(lambda m: str(val.get(int(m.group(1)), 0)), lines[idx].rsplit('=', 1)[0])
    return eval(expr, {'__builtins__': {}})
for idx in [22093, 25539, 2071]:
    print(f"  eq {idx}: value at best = {evline(idx, best)}")

# What does M1=0 need: x_4007 = -L1; with x_5101=1 -> x_30317=-L1
# Print what the three partners must be set to under V=1
print("\n=== Under decoupling x_5101=x_32017=x_26789=1 (V=1): required partners ===")
print(f"  x_30317 = -L1 = {-L1}")
print(f"  x_2936  = 537773*L3 = {537773*L3}")
if L2 % 6672769 == 0:
    print(f"  x_5146  = L2/6672769 = {L2//6672769}")
else:
    print(f"  x_5146  = L2/6672769 NOT integer (need 6672769|L2); L2%6672769={L2%6672769}")
print(f"  and gate outputs: x_4007=x_5101*x_30317, x_29804=x_32017*x_5146, x_35605=x_26789*x_2936")
