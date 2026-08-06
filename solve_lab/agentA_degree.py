#!/usr/bin/env python3
"""Fast mod-p partial evaluator over the S,T cone. Measure the multilinear degree
of S mod p and T mod p as functions of the 78 boolean cone bits."""
import json, sys, random
from agentA_harness import (p, order, gcode, definer, gates, freeinp, anc,
                            backward_cone, load_solution, forward, NVARS)

boolset = set(json.load(open('boolbits.json'))['boolvars'])
S_root, T_root = 35389, 6671

allS, freeS = backward_cone(S_root)
allT, freeT = backward_cone(T_root)
cone = allS | allT
conefree = freeS | freeT
bitsS = sorted(freeS & boolset)
bitsT = sorted(freeT & boolset)
allbits = sorted(conefree & boolset)
slackS = sorted(freeS - boolset)
slackT = sorted(freeT - boolset)
print(f"cone vars={len(cone)} conefree={len(conefree)} allbits={len(allbits)} slack={len(conefree)-len(allbits)}")
print(f"bitsS={len(bitsS)} bitsT={len(bitsT)} overlap={len(set(bitsS)&set(bitsT))}")
print(f"bitsS==bitsT: {bitsS==bitsT}")

# restricted topo order: positions k where order[k] in cone
cone_pos = [k for k, t in enumerate(order) if t in cone]
cone_targets = [order[k] for k in cone_pos]
cgcode = [gcode[k] for k in cone_pos]
print(f"cone gate evals per forward: {len(cone_pos)}")

# baseline
base = load_solution('best/new_instance_partial_39013.json')
forward(base)  # ensure consistent
# reduce cone vars mod p into working array vm (full length, but we only touch cone)
vm = [0]*NVARS
for v in cone:
    vm[v] = base[v] % p
# check current bit values
bitval = {b: base[b] for b in allbits}
nonbin = [b for b in allbits if base[b] not in (0, 1)]
print(f"non-binary bits among cone bits: {len(nonbin)}  {nonbin[:10]}")

ns = {'__builtins__': {}, 'v': vm}
def peval():
    """partial forward-eval mod p over cone; returns (S,T) mod p. vm free inputs must be set."""
    for j, t in enumerate(cone_targets):
        vm[t] = eval(cgcode[j], ns) % p
    return vm[S_root], vm[T_root]

# sanity
S0, T0 = peval()
print(f"peval S0={S0}\n      Sref={base[35389]%p}\n  match={S0==base[35389]%p}")
print(f"peval T0={T0}\n      Tref={base[6671]%p}\n  match={T0==base[6671]%p}")

def set_bits(flips):
    """set cone free inputs to baseline, then flip given bits (0<->1)."""
    for v in conefree:
        vm[v] = base[v] % p
    for b in flips:
        vm[b] = (1 - base[b]) % p  # flip 0<->1

# ---- single-bit differences ----
set_bits([]); Sb, Tb = peval()
assert (Sb, Tb) == (S0, T0)
dS = {}; dT = {}
for b in allbits:
    set_bits([b]); S1, T1 = peval()
    dS[b] = (S1 - S0) % p
    dT[b] = (T1 - T0) % p
nzS = [b for b in allbits if dS[b] != 0]
nzT = [b for b in allbits if dT[b] != 0]
print(f"\nsingle-bit: S sensitive to {len(nzS)}/{len(allbits)} bits; T sensitive to {len(nzT)}/{len(allbits)} bits")

# ---- test degree: second-order differences on all pairs among sensitive bits ----
def second_diff(bi, bj, which):
    set_bits([]); S00, T00 = peval()
    set_bits([bi]); S10, T10 = peval()
    set_bits([bj]); S01, T01 = peval()
    set_bits([bi, bj]); S11, T11 = peval()
    if which == 'S':
        return (S11 - S10 - S01 + S00) % p
    else:
        return (T11 - T10 - T01 + T00) % p

import itertools
sens = sorted(set(nzS) | set(nzT))
print(f"union sensitive bits: {len(sens)}")
quadS = []; quadT = []
pairs = list(itertools.combinations(sens, 2))
print(f"testing {len(pairs)} pairs for 2nd-order term...")
for bi, bj in pairs:
    d = second_diff(bi, bj, 'S')
    if d != 0: quadS.append((bi, bj, d))
    d = second_diff(bi, bj, 'T')
    if d != 0: quadT.append((bi, bj, d))
print(f"nonzero 2nd-order (quadratic) terms: S={len(quadS)}  T={len(quadT)}")
if quadS[:5]: print("  S quad sample:", [(a, b) for a, b, _ in quadS[:8]])
if quadT[:5]: print("  T quad sample:", [(a, b) for a, b, _ in quadT[:8]])

json.dump({'bitsS': bitsS, 'bitsT': bitsT, 'allbits': allbits, 'slackS': slackS, 'slackT': slackT,
           'nzS': nzS, 'nzT': nzT, 'sens': sens,
           'dS': {str(k): v for k, v in dS.items()}, 'dT': {str(k): v for k, v in dT.items()},
           'S0': S0, 'T0': T0,
           'quadS': [[a, b, c] for a, b, c in quadS], 'quadT': [[a, b, c] for a, b, c in quadT]},
          open('agentA_degree.json', 'w'))
print("saved agentA_degree.json")
