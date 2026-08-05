#!/usr/bin/env python3
"""Forward solver v2: zero only the free BOOLEAN inputs; let propagation solve
the value-wires (x_B in huge atoms  bit*(x_B-HUGE)=s*x_C  gets solved once the
bit and x_C are known). Only zero-fill genuinely-free leftovers at the end."""
import json, time
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, substitute, solve_single
from solve_forward import Engine

NVARS = 38748

def boolean_vars(atoms):
    bset = set()
    for poly in atoms:
        vs = atom_vars(poly)
        if len(vs) == 1 and len(poly) == 2:
            v = next(iter(vs))
            if (poly.get((v,)) == 1 and poly.get((v, v)) == -1) or \
               (poly.get((v,)) == -1 and poly.get((v, v)) == 1):
                bset.add(v)
    return bset

def main():
    t0 = time.time()
    atoms = load_atoms()
    eng = Engine(atoms)
    eng.propagate()
    print(f"initial propagation: {eng.n_assigned()} assigned, contra={len(eng.contra)}")

    bset = boolean_vars(atoms)
    print(f"boolean vars: {len(bset)}")

    # Phase 1: zero unassigned booleans one at a time, propagate between
    free_bits = [v for v in bset if eng.val[v] is None]
    print(f"free bits to zero: {len(free_bits)}")
    for v in free_bits:
        if eng.val[v] is None:
            eng.assign(v, 0)
            eng.propagate()
    print(f"after zeroing free bits: {eng.n_assigned()} assigned, contra={len(eng.contra)}")

    # Phase 2: remaining unassigned (value wires not solved, genuinely free) -> 0
    remaining = [v for v in range(NVARS) if eng.val[v] is None]
    print(f"remaining unassigned after bit-zeroing: {len(remaining)}")
    for v in remaining:
        if eng.val[v] is None:
            eng.assign(v, 0)
            eng.propagate()
    print(f"after zero-fill remainder: {eng.n_assigned()} assigned, contra={len(eng.contra)}  time={time.time()-t0:.1f}s")

    v = [x if x is not None else 0 for x in eng.val]
    def ev(poly):
        s = 0
        for m, c in poly.items():
            t = c
            for var in m: t *= v[var]
            s += t
        return s
    violated = [ai for ai, poly in enumerate(atoms) if ev(poly) != 0]
    print(f"violated atoms: {len(violated)}")
    from collections import Counter
    degc = Counter(max((len(m) for m in atoms[ai]), default=0) for ai in violated)
    bigc = sum(1 for ai in violated if any(abs(c) >= 10**20 for c in atoms[ai].values()))
    print(f"violated by degree: {dict(sorted(degc.items()))}, with-big-const: {bigc}")
    json.dump({f"x_{i}": v[i] for i in range(NVARS)}, open('solve_lab/cand_forward2.json', 'w'))
    json.dump(violated, open('solve_lab/violated_forward2.json', 'w'))
    print(f"contradictions ({len(eng.contra)}): {eng.contra[:15]}")
    print("wrote cand_forward2.json")

if __name__ == '__main__':
    main()
