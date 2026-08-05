#!/usr/bin/env python3
"""Cone re-solve with PRODUCT-PRIORITY matching (custom, no SAT/SMT).

The true circuit orients every product/square gate forward: output = the unique
linear variable, inputs = the product operands. best instead used divisions
(x_18274 = x_6773/x_8821), so its whole forward cone is consistent with the wrong
x_9770. We rebuild the gate<->variable matching seeding product/square outputs
FIRST (forced orientation), then augment sum/affine gates, then forward-evaluate
the full ripple cone of the two wrong roots. Exact integers."""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(2_000_000)
NV_PRIM = 4

def build_matching(atoms, avars):
    prim = [len(avars[a]) <= NV_PRIM for a in range(len(atoms))]
    cand_atoms = defaultdict(list); outvars = {}; has_prod = [False]*len(atoms)
    for a in range(len(atoms)):
        if not prim[a]: continue
        bad = set(); lin = set()
        for m in atoms[a]:
            if len(m) == 1: lin.add(m[0])
            else: bad.update(m); has_prod[a] = True
        outvars[a] = lin - bad
        for v in (lin - bad): cand_atoms[v].append(a)
    match = [None]*NVARS; used = [None]*len(atoms)
    # 1) seed product/square gates: their single linear var is the forced output
    for a in range(len(atoms)):
        if prim[a] and has_prod[a] and len(outvars[a]) == 1:
            v = next(iter(outvars[a]))
            if match[v] is None and used[a] is None:
                match[v] = a; used[a] = v
    seeded = sum(1 for m in match if m is not None)
    # 2) augment remaining vars, but never steal a product's forced output
    def augment(v, seen):
        for a in cand_atoms[v]:
            if a in seen: continue
            seen.add(a)
            cur = used[a]
            if cur is None:
                match[v] = a; used[a] = v; return True
            if has_prod[a]:
                continue  # don't displace a product gate's forced output
            if augment(cur, seen):
                match[v] = a; used[a] = v; return True
        return False
    for v in range(NVARS):
        if match[v] is None and cand_atoms[v]:
            augment(v, set())
    print(f"  product-seeded {seeded}, total matched {sum(1 for m in match if m is not None)}", flush=True)
    return match

def main():
    t0 = time.time()
    atoms = load_atoms(); avars = [atom_vars(p) for p in atoms]
    match = build_matching(atoms, avars)

    best = json.load(open('best/best_partial_39013.json'))
    val = [0]*NVARS
    for k, x in best.items():
        i = int(k[2:])
        if i < NVARS: val[i] = x

    succ = defaultdict(list)
    for v in range(NVARS):
        a = match[v]
        if a is None: continue
        for x in avars[a]:
            if x != v: succ[x].append(v)
    roots = [9770, 3183]
    cone = set(roots); dq = deque(roots)
    while dq:
        x = dq.popleft()
        for v in succ[x]:
            if v not in cone: cone.add(v); dq.append(v)
    cone = [v for v in cone if match[v] is not None]
    print(f"cone: {len(cone)} derived vars ({time.time()-t0:.0f}s)", flush=True)

    cset = set(cone); indeg = {v: 0 for v in cone}; adj = defaultdict(list)
    for v in cone:
        for x in avars[match[v]]:
            if x != v and x in cset:
                adj[x].append(v); indeg[v] += 1
    q = deque([v for v in cone if indeg[v] == 0]); topo = []
    while q:
        v = q.popleft(); topo.append(v)
        for u in adj[v]:
            indeg[u] -= 1
            if indeg[u] == 0: q.append(u)
    rest = [v for v in cone if v not in set(topo)]
    order = topo + rest
    print(f"topo {len(topo)}, cyclic {len(rest)}", flush=True)

    def recompute(v):
        a = match[v]; coef = 0; rs = 0
        for m, c in atoms[a].items():
            if m == (v,): coef += c
            else:
                t = c
                for x in m: t *= val[x]
                rs += t
        if coef == 0 or (-rs) % coef != 0: return None
        return (-rs) // coef

    for it in range(60):
        ch = 0
        for v in order:
            r = recompute(v)
            if r is not None and r != val[v]: val[v] = r; ch += 1
        if ch == 0:
            print(f"cone converged after {it+1} sweeps", flush=True); break
    else:
        print("cone did NOT converge in 60 sweeps", flush=True)

    viol = []
    for a in range(len(atoms)):
        s = 0
        for m, c in atoms[a].items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: viol.append(a)
    print(f"TOTAL violated atoms: {len(viol)}  {viol[:30]} ({time.time()-t0:.0f}s)", flush=True)
    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_cone2.json', 'w'))
    print("wrote cand_cone2.json", flush=True)

if __name__ == '__main__':
    main()
