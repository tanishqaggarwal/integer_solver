#!/usr/bin/env python3
"""Try structured/random bit patterns for the 256 core bits, forward-evaluate,
report violated atoms + equation score. Cheap lottery over patterns."""
import json, time
from collections import deque
from propagate import load_atoms, atom_vars
from repair import ProvEngine, boolean_vars

NVARS = 38748

def main():
    atoms = load_atoms()
    bset = boolean_vars(atoms)
    base = ProvEngine(atoms); base.propagate()
    bval = list(base.val); bprov = list(base.prov); bdom = dict(base.domain)
    mainbits = json.load(open('main_comp.json'))['main_bits']

    def run(ones):
        eng = ProvEngine(atoms)
        eng.val = list(bval); eng.prov = list(bprov); eng.domain = dict(bdom)
        eng.wl = deque(); eng.inwl = [False] * len(atoms)
        for b in ones:
            if eng.val[b] is None: eng.assign(b, 1, ('flip', ()))
        eng.propagate()
        for v in [b for b in bset if eng.val[b] is None]:
            if eng.val[v] is None: eng.assign(v, 0, ('free', ())); eng.propagate()
        for v in range(NVARS):
            if eng.val[v] is None: eng.assign(v, 0, ('free', ())); eng.propagate()
        val = [x if x is not None else 0 for x in eng.val]
        viol = 0
        for poly in atoms:
            s = 0
            for m, c in poly.items():
                t = c
                for var in m: t *= val[var]
                s += t
            if s != 0: viol += 1
        return viol, len(eng.contra), val

    tests = []
    tests.append(("all-bits-0 (baseline)", []))
    tests.append(("all-256-bits-1", list(mainbits)))
    res = json.load(open('flip_results.json'))
    improving = [b for (v, b, nc) in res if v < 4]
    tests.append((f"all-{len(improving)}-improving-bits-1", improving))

    # deterministic pseudo-random subsets (index-based, no RNG needed)
    for frac_i, frac in enumerate([2, 3, 4, 8]):
        sub = [b for j, b in enumerate(sorted(mainbits)) if (j * 2654435761) % frac == 0]
        tests.append((f"~1/{frac} subset ({len(sub)} bits)", sub))

    best = None
    for name, ones in tests:
        viol, contra, val = run(ones)
        print(f"{name}: violated_atoms={viol} contra={contra}", flush=True)
        if best is None or viol < best[0]:
            best = (viol, name, val)
    print(f"\nBEST: {best[1]} with {best[0]} violated atoms")
    if best[0] == 0:
        json.dump({f"x_{i}": best[2][i] for i in range(NVARS)}, open('cand_pattern_solved.json', 'w'))
        print("SOLVED -> cand_pattern_solved.json")

if __name__ == '__main__':
    main()
