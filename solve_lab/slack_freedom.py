#!/usr/bin/env python3
"""Are x_26977 and x_9982 (the slacks in a1817/a30378) FREE, or constrained by
other atoms? If free, then a1817/a30378 are trivially satisfiable by choosing the
slack, and the ONLY hard twist constraints are a44271 (x_3183=x_17728) and
a40782 (Q^2=0). That would collapse the obstruction dramatically."""
import json
from collections import defaultdict
from propagate import load_atoms, atom_vars

def fmt(poly, control):
    parts = []
    for m, c in sorted(poly.items(), key=lambda kv: (-len(kv[0]), kv[0])):
        cc = str(c) if abs(c) < 10**9 else f'H{len(str(abs(c)))}'
        mon = '*'.join('x'+str(x)+('#' if x in control else '') for x in m) if m else '1'
        parts.append(f'{cc}*{mon}')
    return ' + '.join(parts)

def main():
    atoms = load_atoms()
    control = set(json.load(open('control_bits.json')))
    best = json.load(open('best/best_partial_39019.json')); bv = {int(k[2:]): v for k, v in best.items()}
    def val(v): return bv.get(v, 0)
    var_atoms = defaultdict(list)
    for a, poly in enumerate(atoms):
        for v in atom_vars(poly):
            var_atoms[v].append(a)

    for t in (26977, 9982, 6773, 17233):
        ats = var_atoms[t]
        print(f"\n===== x_{t}: in {len(ats)} atoms; best={val(t)}; control={t in control} =====")
        for a in ats:
            poly = atoms[a]
            s = 0
            for m, c in poly.items():
                tt = c
                for x in m: tt *= val(x)
                s += tt
            print(f"  a{a} [nv={len(atom_vars(poly))} resid@best={s}]: {fmt(poly,control)[:170]}")

if __name__ == '__main__':
    main()
