#!/usr/bin/env python3
"""Product-priority matching + combo-definers for residue leaves, iterated.

full_forward left residue leaves (defined by combination atoms in best) as fixed
free inputs, so the ripple couldn't close at those boundaries. Here we ALSO let
each such leaf be recomputed by its original combination definer, then iterate
every defined variable to a global fixpoint. Exact integers."""
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
    for a in range(len(atoms)):
        if prim[a] and has_prod[a] and len(outvars[a]) == 1:
            v = next(iter(outvars[a]))
            if match[v] is None and used[a] is None:
                match[v] = a; used[a] = v
    def augment(v, seen):
        for a in cand_atoms[v]:
            if a in seen: continue
            seen.add(a)
            cur = used[a]
            if cur is None:
                match[v] = a; used[a] = v; return True
            if has_prod[a]: continue
            if augment(cur, seen):
                match[v] = a; used[a] = v; return True
        return False
    for v in range(NVARS):
        if match[v] is None and cand_atoms[v]:
            augment(v, set())
    return match

def main():
    t0 = time.time()
    atoms = load_atoms(); avars = [atom_vars(p) for p in atoms]
    match = build_matching(atoms, avars)
    prov = json.load(open('eval_order.json'))['prov']
    # add combo-definers for residue leaves (unmatched vars that best derived from a combo)
    nleaf = 0
    for v in range(NVARS):
        if match[v] is None and v < len(prov) and prov[v] and prov[v][0] >= 0:
            a = prov[v][0]
            if len(avars[a]) > NV_PRIM:
                # combo can define v only if v appears linearly in it
                if any(m == (v,) for m in atoms[a]):
                    match[v] = a; nleaf += 1
    derived = [v for v in range(NVARS) if match[v] is not None]
    print(f"derived {len(derived)} (+{nleaf} leaf-combos), free {NVARS-len(derived)} ({time.time()-t0:.0f}s)", flush=True)

    best = json.load(open('best/best_partial_39013.json'))
    val = [0]*NVARS
    for k, x in best.items():
        i = int(k[2:])
        if i < NVARS: val[i] = x

    order = json.load(open('eval_order.json'))['order']
    dset = set(derived)
    sweep = [v for v in order if v in dset] + [v for v in derived if v not in set(order)]
    gate = {}
    for v in derived:
        a = match[v]; coef = 0; terms = []
        for m, c in atoms[a].items():
            if m == (v,): coef += c
            else: terms.append((c, m))
        gate[v] = (coef, terms)
    def recompute(v):
        coef, terms = gate[v]; rs = 0
        for c, mv in terms:
            t = c
            for x in mv: t *= val[x]
            rs += t
        if coef == 0 or (-rs) % coef != 0: return None
        return (-rs) // coef
    def nviol():
        n = 0
        for a in range(len(atoms)):
            s = 0
            for m, c in atoms[a].items():
                t = c
                for x in m: t *= val[x]
                s += t
            if s: n += 1
        return n
    for rnd in range(30):
        ch = 0
        for v in sweep:
            r = recompute(v)
            if r is not None and r != val[v]: val[v] = r; ch += 1
        nv = nviol() if (ch == 0 or rnd < 2 or rnd % 4 == 0) else '-'
        print(f"round {rnd}: changed={ch} violated={nv} ({time.time()-t0:.0f}s)", flush=True)
        if ch == 0: break
    viol = []
    for a in range(len(atoms)):
        s = 0
        for m, c in atoms[a].items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: viol.append(a)
    print(f"FINAL violated atoms: {len(viol)}  {viol[:40]}", flush=True)
    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_override.json', 'w'))
    print("wrote cand_override.json", flush=True)

if __name__ == '__main__':
    main()
