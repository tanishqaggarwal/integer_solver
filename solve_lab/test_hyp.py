#!/usr/bin/env python3
"""Test fix DIRECTION: change residue leaves x_23268, x_28091 to match best's
x_9770, x_3183 (so gates 27973/27978 hold without moving x_9770/x_3183), then
forward-iterate and count violations. Compare against the other direction."""
import json, sys, time
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(2_000_000)
NV_PRIM = 4

def build_matching(atoms, avars):
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

def main():
    t0 = time.time()
    atoms = load_atoms(); avars = [atom_vars(p) for p in atoms]
    match = build_matching(atoms, avars)
    best = json.load(open('best/best_partial_39013.json'))
    bval = [0]*NVARS
    for k, x in best.items(): bval[int(k[2:])] = x
    order = json.load(open('eval_order.json'))['order']
    derived = [v for v in range(NVARS) if match[v] is not None]
    dset = set(derived)
    sweep = [v for v in order if v in dset] + [v for v in derived if v not in set(order)]
    gate = {}
    for v in derived:
        a = match[v]; coef = 0; terms = []
        for m, c in atoms[a].items():
            if m == (v,): coef += c
            else: terms.append((c, m))
        gate[v] = (coef, terms)

    def run(seed, label):
        val = list(bval)
        for k, x in seed.items(): val[k] = x
        for rnd in range(30):
            ch = 0
            for v in sweep:
                coef, terms = gate[v]; rs = 0
                for c, mv in terms:
                    t = c
                    for x in mv: t *= val[x]
                    rs += t
                if coef and (-rs) % coef == 0:
                    nv = (-rs)//coef
                    if nv != val[v]: val[v] = nv; ch += 1
            if ch == 0: break
        viol = []
        for a in range(len(atoms)):
            s = 0
            for m, c in atoms[a].items():
                t = c
                for x in m: t *= val[x]
                s += t
            if s: viol.append(a)
        print(f"[{label}] violations {len(viol)}: {viol[:20]} ({time.time()-t0:.0f}s)", flush=True)
        return val, viol

    # x_23268 is a free input? check match
    print(f"match[23268]={match[23268]} match[28091]={match[28091]} match[9770]={match[9770]} match[3183]={match[3183]}", flush=True)
    # Direction B: change residue leaves to match best's x_9770/x_3183
    seedB = {23268: bval[9770], 28091: bval[3183]}
    valB, violB = run(seedB, "B: leaves->best x_9770/x_3183")
    if len(violB) == 0:
        json.dump({f"x_{i}": valB[i] for i in range(NVARS)}, open('cand_hypB.json', 'w'))
        print("*** B SOLVES IT ***", flush=True)

if __name__ == '__main__':
    main()
