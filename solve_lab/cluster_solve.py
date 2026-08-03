#!/usr/bin/env python3
"""Solve the bounded conflict cluster by searching its free bits (custom).

The ~240-var cluster around the 4 violated atoms has 30 free inputs; most are
0/1 control bits. Fix the big residue leaves + non-relevant inputs to best,
identify which free bits actually influence the violated atoms, then enumerate
those bits: forward-evaluate the cluster's gate DAG (exact integers, iterated to
fixpoint) and keep the pattern that zeroes every cluster atom."""
import json, sys, time, itertools
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1_000_000)
NV_PRIM = 4
t0 = time.time()

atoms = load_atoms(); avars = [atom_vars(p) for p in atoms]
best = json.load(open('best/best_partial_39013.json'))
bval = [0]*NVARS
for k, x in best.items(): bval[int(k[2:])] = x

# ---- product-priority matching (global) ----
def build_matching():
    prim = [len(avars[a]) <= NV_PRIM for a in range(len(atoms))]
    cand = defaultdict(list); outv = {}; hp = [False]*len(atoms)
    for a in range(len(atoms)):
        if not prim[a]: continue
        bad = set(); lin = set()
        for m in atoms[a]:
            if len(m) == 1: lin.add(m[0])
            else: bad.update(m); hp[a] = True
        outv[a] = lin - bad
        for v in (lin - bad): cand[v].append(a)
    match = [None]*NVARS; used = [None]*len(atoms)
    for a in range(len(atoms)):
        if prim[a] and hp[a] and len(outv[a]) == 1:
            v = next(iter(outv[a]))
            if match[v] is None and used[a] is None: match[v] = a; used[a] = v
    def aug(v, seen):
        for a in cand[v]:
            if a in seen: continue
            seen.add(a)
            if used[a] is None: match[v] = a; used[a] = v; return True
            if hp[a]: continue
            if aug(used[a], seen): match[v] = a; used[a] = v; return True
        return False
    for v in range(NVARS):
        if match[v] is None and cand[v]: aug(v, set())
    return match
match = build_matching()

# ---- cluster = forward cone of the two wrong roots ----
succ = defaultdict(list)
for v in range(NVARS):
    a = match[v]
    if a is None: continue
    for x in avars[a]:
        if x != v: succ[x].append(v)
cluster = set([9770, 3183]); dq = deque([9770, 3183])
while dq:
    x = dq.popleft()
    for v in succ[x]:
        if v not in cluster: cluster.add(v); dq.append(v)
print(f"cluster {len(cluster)} vars ({time.time()-t0:.0f}s)", flush=True)

# derived-in-cluster in topo order + free inputs
derived = [v for v in cluster if match[v] is not None]
free = [v for v in cluster if match[v] is None]
cset = set(cluster)
indeg = {v: 0 for v in derived}; adj = defaultdict(list)
for v in derived:
    for x in avars[match[v]]:
        if x != v and x in set(derived): adj[x].append(v); indeg[v] += 1
q = deque([v for v in derived if indeg[v] == 0]); topo = []
while q:
    v = q.popleft(); topo.append(v)
    for u in adj[v]:
        indeg[u] -= 1
        if indeg[u] == 0: q.append(u)
cyc = [v for v in derived if v not in set(topo)]
order = topo + cyc
print(f"derived {len(derived)} (topo {len(topo)}, cyc {len(cyc)}), free {len(free)}", flush=True)

# atoms to satisfy: every atom that touches the cluster
target_atoms = [a for a in range(len(atoms)) if avars[a] & cset]
print(f"cluster-touching atoms to satisfy: {len(target_atoms)}", flush=True)

gate = {}
for v in derived:
    a = match[v]; coef = 0; terms = []
    for m, c in atoms[a].items():
        if m == (v,): coef += c
        else: terms.append((c, m))
    gate[v] = (coef, terms)

val = list(bval)
def forward():
    for _ in range(len(order) + 5):
        ch = 0
        for v in order:
            coef, terms = gate[v]; rs = 0
            for c, mv in terms:
                t = c
                for x in mv: t *= val[x]
                rs += t
            if coef and (-rs) % coef == 0:
                nv = (-rs) // coef
                if nv != val[v]: val[v] = nv; ch += 1
        if ch == 0: return
def count_viol():
    n = 0
    for a in target_atoms:
        s = 0
        for m, c in atoms[a].items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: n += 1
    return n

# baseline
forward(); base = count_viol()
print(f"baseline cluster violations: {base}", flush=True)

# influence: which free bits (0/1) change the violation set
freebits = [v for v in free if bval[v] in (0, 1)]
bigfree = [v for v in free if bval[v] not in (0, 1)]
print(f"free bits {len(freebits)}, big free (fixed) {len(bigfree)}: {bigfree}", flush=True)
relevant = []
for v in freebits:
    for x in range(NVARS): val[x] = bval[x]
    val[v] = 1 - bval[v]
    forward()
    if count_viol() != base: relevant.append(v)
for x in range(NVARS): val[x] = bval[x]
print(f"relevant free bits ({len(relevant)}): {relevant} ({time.time()-t0:.0f}s)", flush=True)
json.dump(relevant, open('cluster_relevant_bits.json', 'w'))
