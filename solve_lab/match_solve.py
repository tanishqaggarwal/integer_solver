#!/usr/bin/env python3
"""Exact solve via maximum bipartite matching (primitive gates <-> variables).

Each variable is defined by at most one primitive gate; unmatched variables are
free circuit inputs. A maximum matching recovers the straight-line-program
orientation. Free inputs are seeded from the best solution; derived variables are
forward-evaluated in topological order, exactly over Z. Combination atoms never
define anything (they're redundant checks)."""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(2_000_000)

NV_PRIM = 4

def main():
    t0 = time.time()
    atoms = load_atoms()
    avars = [atom_vars(p) for p in atoms]
    prim = [len(avars[a]) <= NV_PRIM for a in range(len(atoms))]
    # candidate output vars per primitive atom + inverse
    cand_atoms = defaultdict(list)   # var -> primitive atoms that can output it
    outvars = {}                     # atom -> set of candidate output vars
    for a in range(len(atoms)):
        if not prim[a]: continue
        bad = set(); lin = set()
        for m in atoms[a]:
            if len(m) == 1: lin.add(m[0])
            else: bad.update(m)
        ov = lin - bad
        outvars[a] = ov
        for v in ov: cand_atoms[v].append(a)
    print(f"{sum(prim)} primitive atoms; {len(cand_atoms)} vars have a primitive definer ({time.time()-t0:.0f}s)", flush=True)

    # seed matching from existing (primitive) provenance
    prov = json.load(open('eval_order.json'))['prov']
    match = [None]*NVARS      # var -> atom
    used = [None]*len(atoms)  # atom -> var
    for v, p in enumerate(prov):
        if p and p[0] >= 0 and prim[p[0]] and v in outvars.get(p[0], ()):
            if used[p[0]] is None:
                match[v] = p[0]; used[p[0]] = v
    seeded = sum(1 for m in match if m is not None)
    print(f"seeded matching: {seeded} vars ({time.time()-t0:.0f}s)", flush=True)

    # Kuhn augment for still-unmatched vars that have a candidate definer
    def augment(v, seen):
        for a in cand_atoms[v]:
            if a in seen: continue
            seen.add(a)
            if used[a] is None or augment(used[a], seen):
                match[v] = a; used[a] = v; return True
        return False
    todo = [v for v in range(NVARS) if match[v] is None and cand_atoms[v]]
    naug = 0
    for v in todo:
        if augment(v, set()): naug += 1
    ndef = sum(1 for m in match if m is not None)
    free = [v for v in range(NVARS) if match[v] is None]
    print(f"augmented {naug}; matched(derived) {ndef}, free inputs {len(free)} ({time.time()-t0:.0f}s)", flush=True)

    # seed free inputs from best
    best = json.load(open('best/best_partial_39013.json'))
    val = [None]*NVARS
    for k, x in best.items():
        i = int(k[2:])
        if i < NVARS: val[i] = x
    bestval = [val[i] if val[i] is not None else 0 for i in range(NVARS)]
    val = [None]*NVARS
    for v in free: val[v] = bestval[v]

    # topological order of derived vars: v depends on other vars of its gate
    derived = [v for v in range(NVARS) if match[v] is not None]
    indeg = {v: 0 for v in derived}
    adj = defaultdict(list)
    dset = set(derived)
    for v in derived:
        for x in avars[match[v]]:
            if x != v and x in dset:
                adj[x].append(v); indeg[v] += 1
    q = deque([v for v in derived if indeg[v] == 0])
    topo = []
    while q:
        v = q.popleft(); topo.append(v)
        for u in adj[v]:
            indeg[u] -= 1
            if indeg[u] == 0: q.append(u)
    cyclic = [v for v in derived if v not in set(topo)]
    print(f"topo: {len(topo)} ordered, {len(cyclic)} cyclic ({time.time()-t0:.0f}s)", flush=True)

    def solve_for(a, v):
        coef = 0; rest = 0
        for m, c in atoms[a].items():
            if m == (v,): coef += c
            else:
                t = c
                for x in m: t *= (val[x] if val[x] is not None else bestval[x])
                rest += t
        if coef == 0 or (-rest) % coef != 0: return None
        return (-rest) // coef

    nd = 0
    for v in topo:
        r = solve_for(match[v], v)
        val[v] = bestval[v] if r is None else r
        if r is None: nd += 1
    for v in cyclic:  # fall back to best for cyclic remainder
        if val[v] is None: val[v] = bestval[v]
    print(f"forward-eval done, {nd} non-divisible ({time.time()-t0:.0f}s)", flush=True)

    viol = []
    for a in range(len(atoms)):
        s = 0
        for m, c in atoms[a].items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: viol.append(a)
    print(f"TOTAL violated atoms: {len(viol)}", flush=True)
    print("  first 40:", viol[:40], flush=True)
    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_match.json', 'w'))
    print(f"wrote cand_match.json ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
