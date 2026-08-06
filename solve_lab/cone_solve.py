#!/usr/bin/env python3
"""Solve the conflict cone with z3. Seed from violated atoms + all wires whose
near-solution value is 'blown up' (>2^300), BFS r hops over the atom<->var
graph, free those vars in z3, fix the rest to cand_forward2. Expand if unsat."""
import json, time, sys
from collections import defaultdict, deque
import z3

NVARS = 38748
RADIUS = int(sys.argv[1]) if len(sys.argv) > 1 else 3
TIMEOUT = int(sys.argv[2]) if len(sys.argv) > 2 else 300

def main():
    t0 = time.time()
    atoms = [json.loads(l)['poly'] for l in open('solve_lab/atoms/poly_atoms.jsonl')]
    atom_vars = []
    var_atoms = defaultdict(list)
    for i, poly in enumerate(atoms):
        vs = set()
        for m, c in poly: vs.update(m)
        atom_vars.append(vs)
        for v in vs: var_atoms[v].append(i)

    cand = json.load(open('solve_lab/cand_forward2.json'))
    fixed = [0] * NVARS
    for k, x in cand.items():
        fixed[int(k[2:])] = x
    viol = json.load(open('solve_lab/violated_forward2.json'))

    # seed vars: violated atoms' vars + all blown-up wires
    seedvars = set()
    for ai in viol:
        seedvars |= atom_vars[ai]
    blown = [v for v in range(NVARS) if abs(fixed[v]) > (1 << 300)]
    seedvars |= set(blown)
    print(f"violated atoms: {len(viol)}, blown-up wires: {len(blown)}, seed vars: {len(seedvars)}")

    # BFS over bipartite graph
    frontier = set(seedvars)
    cone_vars = set(seedvars)
    cone_atoms = set(viol)
    for r in range(RADIUS):
        new_atoms = set()
        for v in frontier:
            for ai in var_atoms[v]:
                if ai not in cone_atoms:
                    new_atoms.add(ai)
        newv = set()
        for ai in new_atoms:
            cone_atoms.add(ai)
            for v in atom_vars[ai]:
                if v not in cone_vars:
                    newv.add(v)
        cone_vars |= newv
        frontier = newv
        print(f"  radius {r+1}: cone_vars={len(cone_vars)}, cone_atoms={len(cone_atoms)}")

    # boolean vars
    bools = set()
    for poly in atoms:
        vs = set()
        for m, c in poly: vs.update(m)
        if len(vs) == 1 and len(poly) == 2:
            v = next(iter(vs))
            d = {tuple(m): c for m, c in poly}
            if (d.get((v,)) == 1 and d.get((v, v)) == -1) or (d.get((v,)) == -1 and d.get((v, v)) == 1):
                bools.add(v)

    zvars = {}
    def zv(v):
        if v not in zvars: zvars[v] = z3.Int(f"x_{v}")
        return zvars[v]

    s = z3.Solver(); s.set("timeout", TIMEOUT * 1000)
    ncons = 0
    for ai in cone_atoms:
        poly = atoms[ai]
        expr = z3.IntVal(0)
        for m, c in poly:
            term = z3.IntVal(c)
            for var in m:
                term = term * (zv(var) if var in cone_vars else z3.IntVal(fixed[var]))
            expr = expr + term
        s.add(expr == 0); ncons += 1
    for v in cone_vars:
        if v in bools:
            s.add(zv(v) >= 0, zv(v) <= 1)
    print(f"z3 model: {len(zvars)} free vars, {ncons} constraints, {sum(1 for v in cone_vars if v in bools)} bits. building {time.time()-t0:.1f}s")
    r = s.check()
    print(f"z3: {r}  ({time.time()-t0:.1f}s)")
    if r == z3.sat:
        model = s.model()
        out = dict(cand)
        for v in cone_vars:
            val = model.eval(zv(v), model_completion=True)
            out[f"x_{v}"] = int(val.as_long())
        json.dump(out, open('solve_lab/cand_cone.json', 'w'))
        print("wrote cand_cone.json")
    else:
        print("cone unsat/unknown at this radius; expand radius")

if __name__ == '__main__':
    main()
