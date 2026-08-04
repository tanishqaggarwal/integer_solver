#!/usr/bin/env python3
"""Pivotal test: does exact-integer forward-eval from an arbitrary free-input assignment
satisfy all 39013 wiring equations over Z (only 20 core fail)? Perturb best and check."""
import json, random, sys
from agentB_setup import load, Env, p, NVARS, override

data = load()
env = Env(data)
best = {int(k[2:]): int(v) for k, v in json.load(open('best/new_instance_partial_39013.json')).items()}

# exact-integer forward + root eval
valz = [0] * NVARS
def forward_int():
    for v in range(NVARS): valz[v] = 0
    for v, x in env.pins.items(): valz[v] = x
    for v in env.freeset: valz[v] = freevals.get(v, 0)
    for v, x in override.items(): valz[v] = x
    for t, pol in env.gate_poly:
        s = 0
        for m, c in pol.items():
            term = c
            for v in m: term *= valz[v]
            s += term
        valz[t] = s
def root_int(i):
    s = 0
    for m, c in env.root_poly[i].items():
        term = c
        for v in m: term *= valz[v]
        s += term
    return s
def count_fail():
    return [i for i in range(len(env.root_poly)) if root_int(i) != 0]

freevals = {v: best[v] for v in env.freeset if v in best}
forward_int()
F0 = count_fail()
print(f"baseline (best): {len(env.root_poly)-len(F0)}/{len(env.root_poly)}  ({len(F0)} fail)")
print("fail idx:", F0[:25])

# perturb random free inputs
random.seed(1)
frees = list(env.freeset)
for npert in (1, 5, 50):
    freevals = {v: best[v] for v in env.freeset if v in best}
    for _ in range(npert):
        h = random.choice(frees)
        freevals[h] = freevals.get(h, 0) + random.randint(1, 1000)
    forward_int()
    F = count_fail()
    print(f"perturb {npert} free inputs: {len(env.root_poly)-len(F)}/{len(env.root_poly)}  ({len(F)} fail)")
    newfail = set(F) - set(F0)
    print(f"    NEW failures (not in baseline 20): {len(newfail)}  {sorted(newfail)[:15]}")
