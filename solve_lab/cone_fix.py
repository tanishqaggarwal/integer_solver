#!/usr/bin/env python3
"""Surgical cone re-solve (custom, no SAT/SMT).

Only x_9770 and x_3183 are wrong (defined by combination atoms in the wrong
order). Their true primitive gates (27973, 27978) force new values, which ripple
through a small forward cone (~94 vars). Recompute ONLY that cone from the
primitive-gate matching; leave all other 38k variables exactly as in the best
solution. Exact integers; verify with the full checker afterward."""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(2_000_000)
NV_PRIM = 4

def build_matching(atoms, avars):
    prim = [len(avars[a]) <= NV_PRIM for a in range(len(atoms))]
    cand_atoms = defaultdict(list); outvars = {}
    for a in range(len(atoms)):
        if not prim[a]: continue
        bad = set(); lin = set()
        for m in atoms[a]:
            if len(m) == 1: lin.add(m[0])
            else: bad.update(m)
        outvars[a] = lin - bad
        for v in (lin - bad): cand_atoms[v].append(a)
    prov = json.load(open('eval_order.json'))['prov']
    match = [None]*NVARS; used = [None]*len(atoms)
    for v, p in enumerate(prov):
        if p and p[0] >= 0 and prim[p[0]] and v in outvars.get(p[0], ()) and used[p[0]] is None:
            match[v] = p[0]; used[p[0]] = v
    def augment(v, seen):
        for a in cand_atoms[v]:
            if a in seen: continue
            seen.add(a)
            if used[a] is None or augment(used[a], seen):
                match[v] = a; used[a] = v; return True
        return False
    for v in [x for x in range(NVARS) if match[x] is None and cand_atoms[x]]:
        augment(v, set())
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
    print(f"topo {len(topo)}, cyclic {len(rest)}: {rest[:20]}", flush=True)

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

    for it in range(40):
        ch = 0
        for v in order:
            r = recompute(v)
            if r is not None and r != val[v]: val[v] = r; ch += 1
        if ch == 0:
            print(f"cone converged after {it+1} sweeps", flush=True); break

    viol = []
    for a in range(len(atoms)):
        s = 0
        for m, c in atoms[a].items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: viol.append(a)
    print(f"TOTAL violated atoms: {len(viol)}  {viol[:30]} ({time.time()-t0:.0f}s)", flush=True)
    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_cone.json', 'w'))
    print("wrote cand_cone.json", flush=True)

if __name__ == '__main__':
    main()
