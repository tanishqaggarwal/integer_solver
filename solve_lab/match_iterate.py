#!/usr/bin/env python3
"""Matching + iterate-to-fixpoint forward evaluation (handles cyclic SCCs).

Same primitive-gate matching as match_solve.py, but instead of a single
topological pass we sweep all derived variables in eval order, recomputing each
from its gate, until values stop changing. The true solution is a simultaneous
fixpoint of every gate; starting from the best (only 2 gates off), the fixes to
x_9770 / x_3183 ripple through their cones and converge. Exact integers."""
import json, time, sys
from collections import defaultdict
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
    atoms = load_atoms()
    avars = [atom_vars(p) for p in atoms]
    match = build_matching(atoms, avars)
    derived = [v for v in range(NVARS) if match[v] is not None]
    free = [v for v in range(NVARS) if match[v] is None]
    print(f"derived {len(derived)}, free {len(free)} ({time.time()-t0:.0f}s)", flush=True)

    best = json.load(open('best/best_partial_39013.json'))
    val = [0]*NVARS
    for k, x in best.items():
        i = int(k[2:])
        if i < NVARS: val[i] = x

    # sweep order: eval order restricted to derived vars (roughly topological)
    order = json.load(open('eval_order.json'))['order']
    dset = set(derived)
    sweep = [v for v in order if v in dset] + [v for v in derived if v not in set(order)]

    # precompile each derived var's gate: (coef_v, list of (coef, othervars))
    gate = {}
    for v in derived:
        a = match[v]; coef = 0; terms = []
        for m, c in atoms[a].items():
            if m == (v,): coef += c
            else: terms.append((c, tuple(x for x in m)))
        gate[v] = (coef, terms)

    def recompute(v):
        coef, terms = gate[v]; rest = 0
        for c, mv in terms:
            t = c
            for x in mv: t *= val[x]
            rest += t
        if coef == 0 or (-rest) % coef != 0: return None
        return (-rest) // coef

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

    for rnd in range(60):
        changed = 0
        for v in sweep:
            r = recompute(v)
            if r is not None and r != val[v]:
                val[v] = r; changed += 1
        nv = nviol() if (changed == 0 or rnd < 3 or rnd % 5 == 0) else None
        print(f"round {rnd}: changed={changed}  violated={nv}  ({time.time()-t0:.0f}s)", flush=True)
        if changed == 0:
            break

    # final exact check
    viol = []
    for a in range(len(atoms)):
        s = 0
        for m, c in atoms[a].items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: viol.append(a)
    print(f"FINAL violated atoms: {len(viol)}  {viol[:30]}", flush=True)
    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_iterate.json', 'w'))
    print(f"wrote cand_iterate.json ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
