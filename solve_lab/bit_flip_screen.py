#!/usr/bin/env python3
"""Single-bit-flip screen (exact propagation, no SAT/SMT).

best sets all free control bits to 0. The true solution flips a few to 1. For
each control bit, set it to 1 (others 0), re-propagate exactly, and count total
violated atoms + the twist atoms. Any bit that reduces the count is a lead;
promising bits are then combined pairwise."""
import json, time, sys
from propagate import load_atoms, atom_vars
from solve_forward import Engine
NVARS = 38748
TWIST = [27973, 27978, 1817, 30378, 44271]

def solve_with(setbits, atoms, bset, control):
    eng = Engine(atoms); eng.propagate()
    for b in setbits:
        if eng.val[b] is None: eng.assign(b, 1); eng.propagate()
    for v in [x for x in control if eng.val[x] is None]:
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    for v in range(NVARS):
        if eng.val[v] is None: eng.assign(v, 0); eng.propagate()
    val = [x if x is not None else 0 for x in eng.val]
    viol = []
    for ai, poly in enumerate(atoms):
        s = 0
        for m, c in poly.items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: viol.append(ai)
    return viol, val

def main():
    t0 = time.time()
    atoms = load_atoms()
    bset = set()
    for poly in atoms:
        vs = atom_vars(poly)
        if len(vs) == 1 and len(poly) == 2:
            v = next(iter(vs))
            if poly.get((v,)) in (1, -1) and poly.get((v, v)) == -poly.get((v,)): bset.add(v)
    control = json.load(open('control_bits.json'))
    base, _ = solve_with([], atoms, bset, control)
    print(f"baseline violated atoms: {len(base)} ({time.time()-t0:.0f}s)", flush=True)
    results = []
    for i, b in enumerate(control):
        viol, _ = solve_with([b], atoms, bset, control)
        nt = sum(1 for a in TWIST if a in viol)
        results.append((len(viol), nt, b))
        if len(viol) < len(base) or nt < 4:
            print(f"  bit {b}: violated {len(viol)} (twist {nt}/5)  <-- improves", flush=True)
        if (i+1) % 40 == 0:
            print(f"  ...{i+1}/{len(control)} tested ({time.time()-t0:.0f}s)", flush=True)
    results.sort()
    print("best 15 single bits (violated, twist, bit):", results[:15], flush=True)
    json.dump(results, open('bit_screen_results.json', 'w'))
    print(f"done ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
