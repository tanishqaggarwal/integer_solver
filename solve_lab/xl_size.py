#!/usr/bin/env python3
"""Measure the ripple subsystem for XL-style linearization. Perturbed vars P = vars
that differ between plain forward-eval and slack-active state. Find all atoms
touching P; count distinct monomials (treat each as an XL variable). If the atom
count and monomial count are both a few thousand, XL (linearize -> GF(P) solve ->
lift) is feasible. Also check CLOSURE: are P's atoms only coupled to P + fixed
vars, or do they pull in the whole system?"""
import json, time
from collections import defaultdict
from confluent_eval5 import build5, make_forward
from slack_active import make_slack_solver
from propagate import atom_vars

def main():
    t0 = time.time()
    A, kind, info, seq0, bestval, ncyc = build5()
    order = json.load(open('eval_order.json'))['order']
    defset = set(v for v in kind if kind[v] != 'const')
    seq = [v for v in order if v in defset and v not in (9770, 3183)]
    seq += [v for v in (9770, 3183) if v in defset]
    seq += [v for v in defset if v not in set(order) and v not in (9770, 3183)]
    solve = make_forward(kind, info, seq, bestval)
    run, seq2 = make_slack_solver(kind, info, seq, bestval)

    v1 = solve(list(bestval), [1858])
    frozen = {24026: v1[18274]-v1[35186], 27116: v1[17728]-v1[1642]}
    v2 = run(list(v1), frozen)
    P = set(v for v in range(len(v1)) if v1[v] != v2[v])
    print(f"perturbed vars |P| = {len(P)}", flush=True)

    var_atoms = defaultdict(list)
    for a, poly in enumerate(A):
        for v in atom_vars(poly): var_atoms[v].append(a)

    # BFS closure: atoms touching P, then vars in those atoms, iterate a few rounds
    frontier = set(P); allvars = set(P); atomset = set()
    for rnd in range(4):
        newatoms = set()
        for v in frontier:
            for a in var_atoms[v]:
                if a not in atomset: newatoms.add(a)
        atomset |= newatoms
        newvars = set()
        for a in newatoms:
            newvars |= atom_vars(A[a])
        frontier = newvars - allvars
        allvars |= newvars
        # count distinct monomials in atomset
        monos = set()
        for a in atomset:
            for m in A[a]: monos.add(m)
        print(f"  round {rnd}: atoms={len(atomset)} vars={len(allvars)} monomials={len(monos)} newvars={len(frontier)} ({time.time()-t0:.0f}s)", flush=True)
        if len(frontier) == 0:
            print("  CLOSED (ripple subsystem is finite/closed)", flush=True); break
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
