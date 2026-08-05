#!/usr/bin/env python3
"""Test single-bit flippability: flip each of the 78 cone bits, full forward-eval,
count how many equations break. Records new (S,T) mod p for each flip."""
import json, sys
from agentA_harness import (p, order, freeinp, backward_cone, load_solution, forward,
                            eval_fails, NEQ, eqvars)

boolset = set(json.load(open('boolbits.json'))['boolvars'])
allS, freeS = backward_cone(35389)
allT, freeT = backward_cone(6671)
allbits = sorted((freeS | freeT) & boolset)

base = load_solution('best/new_instance_partial_39013.json')
forward(base)
F0 = set(eval_fails(base))
print(f"baseline fails: {len(F0)} -> {sorted(F0)}")
S0, T0 = base[35389] % p, base[6671] % p

results = {}
for b in allbits:
    v = base[:]
    v[b] = 1 - base[b]
    forward(v)
    F = set(eval_fails(v))
    broke = F - F0           # newly-failing (were satisfied, now fail)
    fixed = F0 - F           # newly-satisfied (core moved)
    results[b] = {'nfail': len(F), 'broke': len(broke), 'fixed': len(fixed),
                  'S': v[35389] % p, 'T': v[6671] % p,
                  'dS': (v[35389] - base[35389]) % p, 'dT': (v[6671] - base[6671]) % p}
    tag = 'FLIPPABLE' if len(broke) == 0 else ''
    print(f"bit {b:6d} (was {base[b]}): fails={len(F):5d} broke={len(broke):5d} fixed={len(fixed):2d} {tag}")

flippable = [b for b in allbits if results[b]['broke'] == 0]
print(f"\nfreely-flippable single bits (0 non-core breaks): {len(flippable)} -> {flippable}")
json.dump({'flippable': flippable, 'results': {str(k): v for k, v in results.items()},
           'F0': sorted(F0)}, open('agentA_flip.json', 'w'))
print("saved agentA_flip.json")
