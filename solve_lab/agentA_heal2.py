#!/usr/bin/env python3
"""Test Route-A variants that preserve equality-checks, then a smart residual-canceling
greedy healer over free inputs (accept only moves reducing total fails)."""
import json, re, ast, sys
from agentA_harness import (p, order, gcode, definer, gates, freeinp, anc, backward_cone,
                            load_solution, forward, eqcode, eqvars, lines, NEQ, NVARS)
from collections import Counter
sys.setrecursionlimit(1000000)

base = load_solution('best/new_instance_partial_39013.json'); forward(base)
F0 = set(i for i in range(NEQ) if eval(eqcode[i], {'__builtins__': {}, 'v': base}) != 0)

crit = {16742, 14853, 12186}
for r in [35389, 6671]:
    _, fr = backward_cone(r); crit |= fr
QUOT = {30317, 2936, 5146}

def set_quot(v):
    if v[11150] % p == 0: v[30317] = -(v[11150])//p
    if (537773*v[37758]) % p == 0: v[2936] = (537773*v[37758])//p
    if v[25739] % (6672769*p) == 0: v[5146] = v[25739]//(6672769*p)

def count(v):
    ns = {'__builtins__': {}, 'v': v}
    return set(i for i in range(NEQ) if eval(eqcode[i], ns) != 0)

def apply_variant(changes):
    v = base[:]
    for var, val in changes.items(): v[var] = val
    forward(v); set_quot(v)
    return v

print("=== Route A variants (which slack change breaks fewer) ===")
for label, ch in [
    ("A  x16742:=x24908, x14853:=x12186", {16742: base[24908], 14853: base[12186]}),
    ("A2 x16742:=x24908, x12186:=x14853", {16742: base[24908], 12186: base[14853]}),
]:
    v = apply_variant(ch); F = count(v)
    st = (v[35389] % p == 0, v[6671] % p == 0)
    broke = F - F0; fixed = F0 - F
    print(f" {label}: {NEQ-len(F)}/{NEQ} S,T0={st} broke={len(broke)} fixedcore={len(fixed)}")
    print(f"    broke: {sorted(broke)}")

# pick A2 if it breaks fewer, else A
vA = apply_variant({16742: base[24908], 14853: base[12186]}); FA = count(vA)
vA2 = apply_variant({16742: base[24908], 12186: base[14853]}); FA2 = count(vA2)
if len(FA2) <= len(FA):
    v = vA2; base_changes = {16742: base[24908], 12186: base[14853]}
    print("\nusing variant A2")
else:
    v = vA; base_changes = {16742: base[24908], 14853: base[12186]}
    print("\nusing variant A")

# --- smart greedy healer ---
# eq multiplicity for prioritizing private knobs
allvarcount = Counter()
for i in range(NEQ):
    for x in eqvars[i]: allvarcount[x] += 1

VAR = re.compile(r'x_(\d+)')
determined = set(crit) | set(QUOT) | set(base_changes)
def coeff_in_eq(v, i, var):
    old = v[var]; ns = {'__builtins__': {}, 'v': v}
    b = eval(eqcode[i], ns); v[var] = old+1; c = eval(eqcode[i], ns) - b; v[var] = old
    return c

def heal_round(v):
    ns = {'__builtins__': {}, 'v': v}
    F = sorted(i for i in range(NEQ) if eval(eqcode[i], ns) != 0)
    improved = False
    for i in F:
        # residual
        R = eval(eqcode[i], ns)
        if R == 0: continue
        # candidate free inputs in this eq (linear, not determined), prioritize low multiplicity
        cands = sorted((x for x in eqvars[i] if x in freeinp and x not in determined),
                       key=lambda x: allvarcount[x])
        for w in cands[:40]:
            c = coeff_in_eq(v, i, w)
            if c == 0: continue
            if R % c != 0: continue
            delta = -R // c
            old = v[w]; v[w] = old + delta
            forward(v); set_quot(v)
            Fnew = count(v)
            if len(Fnew) < len(F):
                improved = True
                determined.add(w)
                break
            else:
                v[w] = old; forward(v); set_quot(v)
        if improved: break
    return improved

F = count(v)
print(f"start smart heal: {NEQ-len(F)}/{NEQ} ({len(F)} fail)")
for it in range(60):
    if not heal_round(v):
        print(f" no improving move at iter {it}"); break
    F = count(v)
    print(f" iter {it}: {NEQ-len(F)}/{NEQ} ({len(F)} fail)  S,T0={(v[35389]%p==0, v[6671]%p==0)}", flush=True)
    if len(F) == 0: print("SOLVED"); break
F = count(v)
print(f"FINAL: {NEQ-len(F)}/{NEQ} ({len(F)} fail): {sorted(F)[:30]}")
if NEQ-len(F) > 39013:
    out = {f"x_{i}": v[i] for i in range(NVARS) if v[i] != 0}
    json.dump(out, open(f'best_agentA_{NEQ-len(F)}.json', 'w'))
    print(f"SAVED best_agentA_{NEQ-len(F)}.json")
