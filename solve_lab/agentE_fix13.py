#!/usr/bin/env python3
"""For each of the 13 wire=1 unpacking breaks, find rare-partner handles.
Residual R_i = eq value (wire=1). For each free input in the eq's cone, forward-eval coefficient.
Identify absorbable handles (coef divides R_i). Watch collateral on the other 12 + wiring."""
import json, re, sys
from collections import defaultdict
from agentE_harness import Harness, M2MOD, freeinp, wire_set, eqvars, eqcode, lines, gates
from agentE_common import CORE, VAR, p
sys.setrecursionlimit(1000000)

NC13 = [8429, 11166, 11915, 12594, 23869, 25313, 26785, 31400, 32300, 36106, 36767, 37257, 37666]

H = Harness(base_path='best/new_instance_partial_39013.json', wirev=1)
H.forward()
val = H.val; ns = H.ns

def evi(i):
    ns['v'] = val; return eval(eqcode[i], ns)

# residuals
print("=== residuals of the 13 (wire=1) ===")
R = {}
for i in NC13:
    R[i] = evi(i)
    print(f"  eq {i}: R bits={R[i].bit_length()}, sign={'-' if R[i]<0 else '+'}, #vars={len(eqvars[i])}, #free={len(eqvars[i]&freeinp)}")

# For each equation, candidate handles = free inputs in the equation directly (fast), plus their forward effect.
# Compute forward-eval coefficient of each direct free input on THIS equation and count collateral.
print("\n=== per-equation direct free-input handles + forward coeffs ===")
allF = [i for i in range(len(lines)) if evi(i) != 0]
print(f"total failing now: {len(allF)} (core {len([i for i in allF if i in CORE])})")

def coeff_and_collateral(h, targets, delta=1):
    """Perturb free input h by delta, forward, return dict eq->change for targets, and #other eqs changed."""
    old = val[h]
    base = {i: evi(i) for i in targets}
    val[h] = old + delta; H.forward()
    changed = {i: evi(i) - base[i] for i in targets}
    # collateral: sample check on ALL eqs is expensive; count over union of vars? just check full set once
    ns['v'] = val
    newfail = set(i for i in range(len(lines)) if eval(eqcode[i], ns) != 0)
    val[h] = old; H.forward()
    return changed, newfail

base_fail = set(allF)
for i in NC13:
    frees = sorted(eqvars[i] & freeinp)
    # rank handles: those with small footprint
    print(f"\n eq {i} (R bits {R[i].bit_length()}): {len(frees)} direct free inputs")
    good = []
    for h in frees[:60]:
        old = val[h]
        b = evi(i)
        val[h] = old + 1; H.forward(); c = evi(i) - b; val[h] = old; H.forward()
        if c != 0:
            good.append((h, c))
    # show handles whose coeff divides R_i (absorbable alone)
    for h, c in good:
        div = (R[i] % c == 0)
        print(f"    x_{h}: coeff={c} ({'DIVIDES R' if div else 'no'}), R/c={R[i]//c if div else '-'}")
