#!/usr/bin/env python3
"""Print the core cone symbolically: the exact algebra a QUBO must encode."""
import os, sys, json
from collections import deque
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 's9', 'eff'))
import lib as L

LAB = os.path.join(os.path.dirname(__file__), '..')
P = 2**256 - 2**32 - 977


def fmt(atom):
    parts = []
    for mono, c in sorted(L.polys[atom].items(), key=lambda kv: (-len(kv[0]), kv[0])):
        t = '*'.join(f'x{x}' for x in mono) if mono else '1'
        parts.append(f'{c:+d}*{t}' if abs(c) != 1 or not mono else
                     (f'{"+" if c > 0 else "-"}{t}'))
    return ' '.join(parts)


def main():
    state = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    v = L.load(state)
    d = json.load(open(os.path.join(os.path.dirname(__file__), 'core_cone.json')))
    cone_atoms, cone_vars, nz = d['cone_atoms'], d['cone_vars'], d['nonzero_checks']
    free = {t for t in range(L.NVARS) if t not in L.definer}

    print(f"CORE CONE  ({len(cone_vars)} vars, {len(cone_atoms)} atoms)\n")
    print("--- the two failing CHECKS (must be driven to 0) ---")
    for a in nz:
        val = L.evalpoly(L.polys[a], v)
        print(f"  a{a}: {fmt(a)}")
        print(f"        value = {val}   ({val.bit_length()} bits)  "
              f"value/p = {val // P if val % P == 0 else 'NOT divisible by p'}")
    print("\n--- defining GATES in the cone ---")
    for a in sorted(cone_atoms):
        if a in nz:
            continue
        oc = L.atom_out.get(a)
        tag = f"->x{oc[1]}" if oc else "CHECK"
        print(f"  a{a:6d} {tag:>10}: {fmt(a)}")
    print("\n--- cone VARIABLES and their values ---")
    for x in cone_vars:
        kind = 'FREE' if x in free else f'gate a{L.definer[x]}'
        val = v[x]
        s = str(val) if abs(val) < 10**12 else f"{val.bit_length()}-bit"
        extra = ''
        if x in free and val != 0:
            extra = f"   (val mod p = {'0' if val % P == 0 else str(val % P)[:24]+'...'})"
        print(f"  x{x:6d} [{kind:>10}] = {s}{extra}")


if __name__ == '__main__':
    main()
