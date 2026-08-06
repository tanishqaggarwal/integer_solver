#!/usr/bin/env python3
"""Extract exact reduced linear forms (with coefficients) of x_9770,x_18274,x_3183,
x_17728, and trace the key free wires {31434,34236,35846,26977,27912,28035,6236,10466}:
their defining atom, kind, and best_partial value. The whole twist obstruction lives
in these ~8 wires."""
import json, time
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
from twist_forms import build_pivots, reduce_var

def main():
    t0 = time.time()
    atoms = load_atoms()
    best = json.load(open('best/best_partial_39019.json'))
    bv = [0]*NVARS
    for k, x in best.items(): bv[int(k[2:])] = x
    pivot_of = build_pivots(atoms)
    import sys; sys.setrecursionlimit(100000)
    memo = {}
    P = (1 << 61) - 1
    for w in (9770, 18274, 3183, 17728):
        f, c = reduce_var(w, pivot_of, memo)
        # print coefficients (as signed, mod P -> map to small if near 0/P)
        def sgn(x):
            return x - P if x > P//2 else x
        terms = sorted(f.items())
        print(f"\nx_{w} = " + " + ".join(f"({sgn(cf)})*x_{v}" for v, cf in terms) + f"  [+{sgn(c)}]")
        print(f"   check value: x_{w} best={bv[w]}")

    # trace the key wires
    var_def = {}
    for a, poly in enumerate(atoms):
        pass
    keys = [31434, 34236, 35846, 26977, 27912, 28035, 6236, 10466]
    # find atoms where each key var appears, identify a plausible 'definition'
    appear = defaultdict(list)
    for a, poly in enumerate(atoms):
        for v in atom_vars(poly):
            if v in keys: appear[v].append(a)
    for v in keys:
        val = bv[v]
        vs = '0' if val == 0 else (f'small({val})' if abs(val) < 10**8 else f'HUGE({len(str(abs(val)))}d)')
        print(f"\nx_{v} = {vs}; in {len(appear[v])} atoms: {appear[v][:8]}")
        # show the shortest atom containing it (likely its definition)
        shortest = min(appear[v], key=lambda a: len(atoms[a]))
        poly = atoms[shortest]
        tt = []
        for m, c in poly.items():
            cc = c if abs(c) < 10**12 else f"HUGE({len(str(abs(c)))}d)"
            tt.append(f"{cc}*{'*'.join('x_'+str(x) for x in m) if m else '1'}")
        print(f"   shortest atom {shortest}: {' + '.join(tt)[:160]}")
    print(f"\ndone ({time.time()-t0:.0f}s)", flush=True)

if __name__ == '__main__':
    main()
