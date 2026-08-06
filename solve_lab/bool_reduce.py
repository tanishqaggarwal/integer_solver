#!/usr/bin/env python3
"""Reduce atoms modulo x^2=x for boolean variables, then propagate. This
linearizes boolean sub-expressions and can force many bits that ordinary
propagation misses (e.g. x*(x-1)+4*(x-1) collapses to 4(x-1) -> x=1)."""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, substitute, solve_single
from repair import ProvEngine, boolean_vars

NVARS = 38748

def reduce_bool(poly, bset):
    """Collapse repeated boolean vars in each monomial (x^k -> x)."""
    out = defaultdict(int)
    for m, c in poly.items():
        newm = []
        seen = set()
        for v in sorted(m):
            if v in bset:
                if v in seen:
                    continue
                seen.add(v)
                newm.append(v)
            else:
                newm.append(v)
        out[tuple(sorted(newm))] += c
    return {m: c for m, c in out.items() if c != 0}

def main():
    t0 = time.time()
    atoms = load_atoms()
    bset = boolean_vars(atoms)
    print(f"{len(atoms)} atoms, {len(bset)} boolean vars", file=sys.stderr)
    # reduced atoms
    ratoms = [reduce_bool(p, bset) for p in atoms]
    # drop trivial (empty) atoms
    ratoms = [p for p in ratoms if p]
    print(f"reduced atoms: {len(ratoms)} (from {len(atoms)}) in {time.time()-t0:.1f}s", file=sys.stderr)

    # propagate on reduced atoms; also add boolean domain constraint by keeping x*(x-1) info:
    # we re-add explicit boolean atoms so domains are known
    for b in bset:
        ratoms.append({(b, b): 1, (b,): -1})  # keep boolean constraint (won't reduce here)
    # But reduce_bool would collapse these; add AFTER reduction so they stay as x^2-x.

    eng = ProvEngine(ratoms)
    eng.propagate()
    mainbits = set(json.load(open('main_comp.json'))['main_bits'])
    fb = [b for b in mainbits if eng.val[b] is not None]
    print(f"after boolean-reduced propagation: {eng.n_assigned()} assigned, contra {len(eng.contra)}")
    print(f"MAIN BITS FORCED: {len(fb)} / 256")
    from collections import Counter
    print("  forced main-bit values:", dict(Counter(eng.val[b] for b in fb)))
    allb_forced = [b for b in bset if eng.val[b] is not None]
    print(f"ALL bits forced: {len(allb_forced)} / {len(bset)}")

    # now complete: use these forced bits, propagate ORIGINAL atoms, zero-fill rest, check
    eng2 = ProvEngine(atoms)
    # seed with forced values from reduced propagation
    forced = {i: eng.val[i] for i in range(NVARS) if eng.val[i] is not None}
    eng2.propagate()
    for i, v in forced.items():
        if eng2.val[i] is None:
            eng2.assign(i, v, ('seed', ()))
    eng2.propagate()
    for b in [x for x in bset if eng2.val[x] is None]:
        if eng2.val[b] is None: eng2.assign(b, 0, ('f', ())); eng2.propagate()
    for v in range(NVARS):
        if eng2.val[v] is None: eng2.assign(v, 0, ('f', ())); eng2.propagate()
    val = [x if x is not None else 0 for x in eng2.val]
    viol = 0
    for poly in atoms:
        s = 0
        for m, c in poly.items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s != 0: viol += 1
    print(f"\nAfter seeding forced bits + completing: {viol} violated atoms, {len(eng2.contra)} contra")
    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_boolreduce.json', 'w'))
    json.dump({str(b): eng.val[b] for b in fb}, open('forced_main_bits.json', 'w'))
    print("wrote cand_boolreduce.json, forced_main_bits.json")

if __name__ == '__main__':
    main()
