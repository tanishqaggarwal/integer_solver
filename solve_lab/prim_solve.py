#!/usr/bin/env python3
"""Primitives-only exact re-derivation.

Key fix: only PRIMITIVE gate atoms (<=4 vars) may DEFINE a variable. Combination
atoms (linear combos of gate residuals) are redundant checks and must never
assign. Free inputs (vars no primitive gate can output) are seeded from the best
solution; every derivable variable is forward-computed from its primitive gate.
Exact integer arithmetic throughout."""
import json, time
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, substitute, solve_single, NVARS

NV_PRIM = 4

def main():
    t0 = time.time()
    atoms = load_atoms()
    avars = [atom_vars(p) for p in atoms]
    prim = [len(avars[a]) <= NV_PRIM for a in range(len(atoms))]
    print(f"loaded {len(atoms)} atoms; {sum(prim)} primitive, {len(atoms)-sum(prim)} combination ({time.time()-t0:.0f}s)", flush=True)

    # candidate-output vars: appear only in degree-1 monomials of a primitive atom
    derivable = set()
    for a in range(len(atoms)):
        if not prim[a]: continue
        bad = set(); lin = set()
        for m in atoms[a]:
            if len(m) == 1: lin.add(m[0])
            else: bad.update(m)
        derivable |= (lin - bad)
    free_inputs = [v for v in range(NVARS) if v not in derivable]
    print(f"derivable vars: {len(derivable)}, free inputs: {len(free_inputs)}", flush=True)

    # seed free inputs from best solution
    best = json.load(open('best/best_partial_39013.json'))
    bestval = [0]*NVARS
    for k, x in best.items(): bestval[int(k[2:])] = x
    val = [None]*NVARS
    for v in free_inputs: val[v] = bestval[v]

    # propagate: only primitive atoms may assign
    var_atoms = defaultdict(list)
    for a in range(len(atoms)):
        if prim[a]:
            for v in avars[a]: var_atoms[v].append(a)
    wl = deque(a for a in range(len(atoms)) if prim[a])
    inwl = [prim[a] for a in range(len(atoms))]
    contra = []
    def assign(v, x):
        if val[v] is not None:
            if val[v] != x: contra.append((v, val[v], x))
            return
        val[v] = x
        for a in var_atoms[v]:
            if not inwl[a]: inwl[a] = True; wl.append(a)
    steps = 0
    while wl:
        a = wl.popleft(); inwl[a] = False; steps += 1
        poly = substitute(atoms[a], val)
        uv = atom_vars(poly)
        if len(uv) == 1:
            kind, data = solve_single(poly)
            if kind == 'val': assign(*data)
            elif kind == 'dom':
                u, roots = data
                if len(roots) == 1: assign(u, next(iter(roots)))
        # len 0 or >=2: skip (checks handled at the end)
    ndet = sum(1 for x in val if x is not None)
    print(f"primitive fixpoint: {ndet}/{NVARS} determined, steps={steps}, contradictions={len(contra)} ({time.time()-t0:.0f}s)", flush=True)
    for c in contra[:8]: print("   contra", str(c)[:100], flush=True)

    # fill any still-undetermined from best (residual free vars)
    nfill = 0
    for v in range(NVARS):
        if val[v] is None: val[v] = bestval[v]; nfill += 1
    print(f"filled {nfill} remaining undetermined from best", flush=True)

    # full exact violation check
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
    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_prim.json', 'w'))
    print(f"wrote cand_prim.json ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
