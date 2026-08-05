#!/usr/bin/env python3
"""Minimal augmenting-chain fix over best's prov orientation (no SAT/SMT).

39019 satisfies the 2 sum-gates but violates the compensation checks 1817
(x_18274=x_9770) and 44271/30378 (x_17728=x_3183). Fix: re-orient the minimal
chain of best's prov definers so these checks become definers, cascading
1817->x_18274, 4954->x_6773, 13202->x_26517, 13203->x_15690, dropping the
terminal combination atom. Recompute the chain vars exactly; iterate for cycles."""
import json, sys, time
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1_000_000)
NV_PRIM = 4

def main():
    t0 = time.time()
    atoms = load_atoms(); avars = [atom_vars(p) for p in atoms]
    prim = [len(avars[a]) <= NV_PRIM for a in range(len(atoms))]
    outv = {}
    for a in range(len(atoms)):
        bad = set(); lin = set()
        for m in atoms[a]:
            if len(m) == 1: lin.add(m[0])
            else: bad.update(m)
        outv[a] = lin - bad
    cand_atoms = defaultdict(list)
    for a in range(len(atoms)):
        if prim[a]:
            for v in outv[a]: cand_atoms[v].append(a)

    best = json.load(open('best/best_partial_39019.json'))
    val = [0]*NVARS
    for k, x in best.items(): val[int(k[2:])] = x
    prov = json.load(open('eval_order.json'))['prov']
    df = [None]*NVARS; owner = [None]*len(atoms)
    for v, p in enumerate(prov):
        if p and p[0] >= 0: df[v] = p[0]; owner[p[0]] = v

    changed_def = set()
    LOCKED = {9770, 3183}   # defined by their true sum-gates 27973/27978 — never displace
    def augment(atom):
        """atom must claim a linear output var via alternating path to combo/free."""
        def dfs(a, seen):
            for v in outv[a]:
                if v in seen or v in LOCKED: continue
                seen.add(v)
                cur = df[v]
                if cur == a: continue
                if cur is None:
                    df[v] = a; owner[a] = v; changed_def.add(v); return True
                if not prim[cur]:                       # drop redundant combo
                    owner[cur] = None; df[v] = a; owner[a] = v; changed_def.add(v); return True
                if dfs(cur, seen):
                    df[v] = a; owner[a] = v; changed_def.add(v); return True
            return False
        return dfs(atom, set())

    # the violated primitive checks in 39019
    checks = [1817, 30378, 44271]
    for a in checks:
        if owner[a] is None:
            ok = augment(a)
            print(f"augment {a}: {'ok' if ok else 'FAIL'} -> defines x_{owner[a]}", flush=True)

    print(f"changed definitions: {sorted(changed_def)} ({time.time()-t0:.0f}s)", flush=True)

    # recompute changed-def vars (and any var whose def changed) in topo order, iterate for cycles
    recompute_set = set(changed_def)
    order_list = list(recompute_set)
    def recompute(v):
        a = df[v]; coef = 0; rs = 0
        for m, c in atoms[a].items():
            if m == (v,): coef += c
            else:
                t = c
                for x in m: t *= val[x]
                rs += t
        if coef == 0 or (-rs) % coef != 0: return None
        return (-rs) // coef
    for _ in range(200):
        ch = 0
        for v in order_list:
            r = recompute(v)
            if r is not None and r != val[v]: val[v] = r; ch += 1
        if ch == 0: break

    viol = []
    for a in range(len(atoms)):
        s = 0
        for m, c in atoms[a].items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: viol.append(a)
    print(f"violated atoms: {len(viol)}: {viol[:30]} ({time.time()-t0:.0f}s)", flush=True)
    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_chainfix.json', 'w'))
    print("wrote cand_chainfix.json", flush=True)

if __name__ == '__main__':
    main()
