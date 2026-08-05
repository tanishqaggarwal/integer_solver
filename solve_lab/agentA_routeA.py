#!/usr/bin/env python3
"""Route A: set x_16742:=x_24908 (x_3558=0) and x_14853:=x_12186 (x_29322=0) EXACTLY,
forward-eval, set quotient handles, measure damage. Then attempt heal."""
import json, sys
from agentA_harness import (p, load_solution, forward, eval_fails, NEQ, order,
                            gates, definer, freeinp, eqvars)

base = load_solution('best/new_instance_partial_39013.json'); forward(base)
F0 = set(eval_fails(base))
print(f"baseline: {NEQ-len(F0)}/{NEQ} ({len(F0)} fail)")

def try_variant(changes, label):
    v = base[:]
    for var, val in changes.items():
        v[var] = val
    forward(v)
    # set quotient handles
    for hv, expr in [(30317, lambda: -(v[11150])//p if v[11150] % p == 0 else None),
                     (2936, lambda: (537773*v[37758])//p if (537773*v[37758]) % p == 0 else None),
                     (5146, lambda: v[25739]//(6672769*p) if v[25739] % (6672769*p) == 0 else None)]:
        r = expr()
        if r is not None: v[hv] = r
    F = set(eval_fails(v))
    Snew = v[35389] % p; Tnew = v[6671] % p
    print(f"\n[{label}] S mod p={'0' if Snew==0 else 'NONZERO'}  T mod p={'0' if Tnew==0 else 'NONZERO'}")
    print(f"  x_3558={v[3558]} (mod p={v[3558]%p})  x_29322={v[29322]} (mod p={v[29322]%p})")
    print(f"  x_11150 % p = {v[11150]%p==0}  x_37758*537773 % p = {(537773*v[37758])%p==0}  x_25739 % (6672769p) = {v[25739]%(6672769*p)==0}")
    print(f"  satisfied {NEQ-len(F)}/{NEQ} ({len(F)} fail)")
    broke = F - F0; fixed = F0 - F
    print(f"  broke {len(broke)} (were ok): {sorted(broke)[:25]}")
    print(f"  fixed {len(fixed)} (core): {sorted(fixed)}")
    return v, F

# Route A: x_3558=0 exactly, x_29322=0 exactly
v, F = try_variant({16742: base[24908], 14853: base[12186]}, "Route A (x16742:=x24908, x14853:=x12186)")

# also try: only x_3558=0 (leave x_29322 as is), general S,T
# with x_29322 nonzero and x_3558=0: S = a*d^2 - 0 = a*d^2 (nonzero), so S!=0. skip.

# Save if improved
if len(F) < len(F0):
    out = {f"x_{i}": v[i] for i in range(len(v)) if v[i] != 0}
    json.dump(out, open('agentA_routeA_result.json', 'w'))
    print("saved agentA_routeA_result.json")
