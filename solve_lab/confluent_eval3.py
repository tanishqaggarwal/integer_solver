#!/usr/bin/env python3
"""Confluent forward evaluator with residue loading AND correct topo order.

Fixes the ordering bug: best's eval_order put x_9770 before its sum-gate input
x_35186. We rebuild a proper topological order for the overridden orientation
(x_9770<-27973, x_3183<-27978, huge-atom targets load conditionally on their bit).
Then a bit flip's residue load propagates correctly to the check gates. Screen
singles + relevant pairs; verify hits in Z."""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
from modp import P, inv
sys.setrecursionlimit(2_000_000)

def boolean_vars(atoms):
    bset = set()
    for poly in atoms:
        if len(poly) == 2:
            ms = sorted(poly.keys(), key=len)
            if len(ms[0]) == 1 and len(ms[1]) == 2 and ms[1] == (ms[0][0], ms[0][0]) and poly[ms[0]] == -poly[ms[1]]:
                bset.add(ms[0][0])
    return bset

def build():
    atoms = load_atoms(); avars = [atom_vars(p) for p in atoms]
    bset = boolean_vars(atoms)
    load = {}
    for a, poly in enumerate(atoms):
        prod2 = [(m, c) for m, c in poly.items() if len(m) == 2 and m[0] != m[1]]
        if len(prod2) != 1: continue
        (m2, c2) = prod2[0]
        bigcoef = any(abs(c) > 10**40 for m, c in poly.items() if len(m) == 1)
        for bit, xB in [(m2[0], m2[1]), (m2[1], m2[0])]:
            if bit in bset and bigcoef and xB not in load:
                load[xB] = (a, bit, c2 % P); break
    prov = json.load(open('eval_order.json'))['prov']
    best = json.load(open('best/best_partial_39019.json'))
    bestval = [0]*NVARS
    for k, x in best.items(): bestval[int(k[2:])] = x
    df = [None]*NVARS
    for v, p in enumerate(prov):
        if p and p[0] >= 0: df[v] = p[0]
    df[9770] = 27973; df[3183] = 27978
    for xB in load: df[xB] = load[xB][0]

    gate = {}; deps = {}
    for v in range(NVARS):
        if v in load:
            a, bit, cbx = load[v]
            terms = [(c % P, m) for m, c in atoms[a].items() if not (len(m) == 2 and v in m and bit in m)]
            gate[v] = ('load', bit, cbx, terms)
            d = set([bit])
            for c, m in terms: d.update(m)
            d.discard(v); deps[v] = d
            continue
        a = df[v]
        if a is None: continue
        coef = 0; terms = []
        for m, c in atoms[a].items():
            if m == (v,): coef = (coef + c) % P
            else: terms.append((c % P, m))
        gate[v] = ('gate', coef, terms)
        d = set(avars[a]); d.discard(v); deps[v] = d

    # topological sort (Kahn) over derived vars
    derived = set(gate)
    indeg = {v: 0 for v in derived}; adj = defaultdict(list)
    for v in derived:
        for x in deps[v]:
            if x in derived:
                adj[x].append(v); indeg[v] += 1
    q = deque([v for v in derived if indeg[v] == 0]); topo = []
    while q:
        v = q.popleft(); topo.append(v)
        for u in adj[v]:
            indeg[u] -= 1
            if indeg[u] == 0: q.append(u)
    cyc = [v for v in derived if v not in set(topo)]
    seq = topo + cyc
    return atoms, gate, seq, [x % P for x in bestval], len(cyc)

def main():
    t0 = time.time()
    atoms, gate, seq, bm, ncyc = build()
    print(f"built: {len(seq)} derived, {ncyc} cyclic ({time.time()-t0:.0f}s)", flush=True)
    CHECKS = [1817, 30378, 44271]
    cp = [[(c % P, m) for m, c in atoms[a].items()] for a in CHECKS]

    def forward(setbits, extra_rounds=0):
        val = list(bm)
        for b in setbits: val[b] = 1
        for _ in range(1 + extra_rounds):
            for v in seq:
                g = gate[v]
                if g[0] == 'load':
                    _, bit, cbx, terms = g
                    if val[bit] == 0: val[v] = 0
                    else:
                        rs = 0
                        for c, m in terms:
                            t = c
                            for x in m: t = (t * val[x]) % P
                            rs = (rs + t) % P
                        val[v] = (-rs * inv((cbx * val[bit]) % P)) % P
                else:
                    _, coef, terms = g
                    if coef == 0: continue
                    rs = 0
                    for c, m in terms:
                        t = c
                        for x in m: t = (t * val[x]) % P
                        rs = (rs + t) % P
                    val[v] = (-rs * inv(coef)) % P
        res = []
        for poly in cp:
            s = 0
            for c, m in poly:
                t = c
                for x in m: t = (t * val[x]) % P
                s = (s + t) % P
            res.append(s)
        return res, val

    base, _ = forward([])
    print(f"baseline checks (0=pass): {[1 if x else 0 for x in base]} ({time.time()-t0:.0f}s)", flush=True)

    control = json.load(open('control_bits.json'))
    changed = []
    for b in control:
        r, _ = forward([b])
        if r != base:
            npass = sum(1 for x in r if x == 0)
            changed.append((npass, b, [1 if x else 0 for x in r]))
    changed.sort(reverse=True)
    print(f"single bits that CHANGE checks: {len(changed)} ({time.time()-t0:.0f}s)", flush=True)
    for npass, b, r in changed[:30]:
        print(f"  bit {b}: pass {npass}/3  {r}", flush=True)
    json.dump([[n, b, r] for n, b, r in changed], open('confluent3_single.json', 'w'))

if __name__ == '__main__':
    main()
