#!/usr/bin/env python3
"""Check the equality-chain hypothesis: checks want x_24908=x_16742=x_19083 and
x_12186=x_14853=x_1308. Inspect partner gates x_19083, x_1308 and their drivers."""
import json
from agentA_harness import (p, load_solution, forward, gates, definer, order, freeinp,
                            backward_cone)

base = load_solution('best/new_instance_partial_39013.json'); forward(base)
gdef = {t: gates[definer[t]][1] for t in order}

def info(v):
    allv, fr = backward_cone(v)
    return gdef.get(v, '(free)'), len(fr), sorted(fr)[:8]

print("=== baseline values (mod p) ===")
for v in [24908, 16742, 19083, 12186, 14853, 1308]:
    print(f"  x_{v} = ...{base[v]%p % 10**20}  (free={v in freeinp})")

print("\n=== equality checks at baseline ===")
print(f"  x_24908 - x_16742 (=x_3558)  mod p = {(base[24908]-base[16742])%p}")
print(f"  x_16742 - x_19083            mod p = {(base[16742]-base[19083])%p}")
print(f"  x_24908 - x_19083            mod p = {(base[24908]-base[19083])%p}")
print(f"  x_14853 - x_12186 (=x_29322) mod p = {(base[14853]-base[12186])%p}")
print(f"  x_14853 - x_1308             mod p = {(base[14853]-base[1308])%p}")
print(f"  x_12186 - x_1308             mod p = {(base[12186]-base[1308])%p}")
print(f"  exact: x_16742==x_19083? {base[16742]==base[19083]}   x_14853==x_1308? {base[14853]==base[1308]}")
print(f"  exact: x_24908==x_19083? {base[24908]==base[19083]}   x_12186==x_1308? {base[12186]==base[1308]}")

print("\n=== partner gate definitions & cones ===")
for v in [19083, 1308, 24908, 12186]:
    rhs, nfr, sample = info(v)
    print(f"  x_{v}: rhs={rhs[:80]}")
    print(f"        free-cone size={nfr} sample={sample}")

# trace x_19083 driver chain
print("\n=== trace x_19083 backward ===")
def trace(v, depth=0, seen=None):
    if seen is None: seen=set()
    if v in seen or depth>6: return
    seen.add(v)
    r = gdef.get(v, '(FREE)')
    vids = gates[definer[v]][2] if v in definer else ()
    print("  "*depth + f"x_{v} = {r[:60]}  vids={list(vids)[:6]}")
    if v in freeinp: return
    for u in vids[:3]:
        trace(u, depth+1, seen)
trace(19083)
print("--- x_1308 ---")
trace(1308)
