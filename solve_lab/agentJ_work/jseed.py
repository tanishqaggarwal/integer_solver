#!/usr/bin/env python3
"""Greedy seeding: how few free inputs suffice to determine everything?"""
import os, pickle, time, heapq
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
P = pickle.load(open(os.path.join(HERE, 'jpoly.pkl'), 'rb'))
polys, defcands = P['polys'], P['defcands']
NA = len(polys)
NV = 38748

varsof = []
for p in polys:
    s = set()
    for k in p:
        s.update(k)
    varsof.append(s)
occ = defaultdict(list)
for i, s in enumerate(varsof):
    for v in s:
        occ[v].append(i)

known = set()
unk = [len(varsof[i]) for i in range(NA)]
fired = [False] * NA
order = []
q = []

def push_ready():
    for i in range(NA):
        if not fired[i] and unk[i] == 1:
            q.append(i)

def run():
    while q:
        i = q.pop()
        if fired[i] or unk[i] != 1:
            continue
        miss = varsof[i] - known
        if len(miss) != 1:
            continue
        v = next(iter(miss))
        if v not in defcands[i]:
            continue
        fired[i] = True
        order.append(('def', i, v))
        known.add(v)
        for j in occ[v]:
            if not fired[j]:
                unk[j] -= 1
                if unk[j] == 1:
                    q.append(j)

def declare(v, tag):
    if v in known:
        return
    known.add(v)
    order.append((tag, None, v))
    for j in occ[v]:
        if not fired[j]:
            unk[j] -= 1
            if unk[j] == 1:
                q.append(j)
    run()

push_ready(); run()
print("after pins:", len(known), "known;", sum(fired), "fired")

# booleans: atoms v^2 - v  (single variable, quadratic)
bools = []
for i in range(NA):
    if len(varsof[i]) == 1:
        p = polys[i]
        if any(len(k) == 2 for k in p):
            bools.append((i, next(iter(varsof[i]))))
print("single-var quadratic atoms:", len(bools))
for i, v in bools:
    if v not in known:
        declare(v, 'bool')
print("after booleans:", len(known), "known;", sum(fired), "fired")

# greedy: repeatedly declare the unknown variable that unlocks the most
t0 = time.time()
seeds = []
while True:
    remaining = [i for i in range(NA) if not fired[i] and len(varsof[i] - known) >= 1]
    if not remaining:
        break
    # candidate seeds: unknown vars in atoms with 2 unknowns, ranked by degree
    cnt = Counter()
    for i in remaining:
        m = varsof[i] - known
        if len(m) == 2:
            for v in m:
                cnt[v] += 1
    if not cnt:
        for i in remaining:
            for v in varsof[i] - known:
                cnt[v] += 1
    if not cnt:
        break
    v = max(cnt, key=lambda z: (cnt[z], len(occ[z])))
    seeds.append(v)
    declare(v, 'seed')
    if len(seeds) % 200 == 0:
        print(f"  seeds={len(seeds)} known={len(known)} fired={sum(fired)} ({time.time()-t0:.0f}s)")

print("SEEDS NEEDED:", len(seeds))
print("known:", len(known), "fired:", sum(fired), "unfired:", NA - sum(fired))
unfired = [i for i in range(NA) if not fired[i]]
print("unfired with 0 unknowns:", sum(1 for i in unfired if not (varsof[i] - known)))
pickle.dump({'order': order, 'seeds': seeds, 'known': known, 'fired': fired},
            open(os.path.join(HERE, 'jseed.pkl'), 'wb'))
