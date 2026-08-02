#!/usr/bin/env python3
"""Isolate the giant residual component and characterize how its free bits are
constrained: dump the huge-constant atoms and the atoms linear in bits."""
import json
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, substitute

NVARS = 38748

def main():
    atoms = load_atoms()
    part = json.load(open('solve_lab/partial_assignment.json'))
    val = [None] * NVARS
    for k, x in part.items():
        val[int(k[2:])] = x

    bool_vars = set()
    for poly in atoms:
        d = {tuple(m): c for m, c in poly.items()}
        vs = atom_vars(poly)
        if len(vs) == 1 and len(d) == 2:
            v = next(iter(vs))
            if (d.get((v,)) == 1 and d.get((v, v)) == -1) or (d.get((v,)) == -1 and d.get((v, v)) == 1):
                bool_vars.add(v)

    # residual atoms and union-find
    unsolved = []
    for ai, poly in enumerate(atoms):
        red = substitute(poly, val)
        uv = atom_vars(red)
        if len(uv) >= 2:
            unsolved.append((ai, red, uv))
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for ai, red, uv in unsolved:
        uv = list(uv)
        for v in uv: find(v)
        for i in range(1, len(uv)): union(uv[0], uv[i])
    comp = defaultdict(set)
    for v in list(parent.keys()):
        comp[find(v)].add(v)
    main_root = max(comp, key=lambda r: len(comp[r]))
    main_vars = comp[main_root]
    main_bits = sorted(v for v in main_vars if v in bool_vars)
    print(f"main component: {len(main_vars)} vars, {len(main_bits)} bits")

    # atoms fully inside main component
    main_atoms = [(ai, red, uv) for ai, red, uv in unsolved if uv & main_vars]
    print(f"main-component residual atoms: {len(main_atoms)}")

    # classify: huge-const atoms, degree, bit-linear
    huge = []
    bitlin = []   # atoms linear whose vars are all bits
    deg_hist = defaultdict(int)
    for ai, red, uv in main_atoms:
        deg = max((len(m) for m in red), default=0)
        deg_hist[deg] += 1
        if any(abs(c) >= 10**20 for c in red.values()):
            huge.append((ai, red, uv))
        if deg == 1 and uv <= set(main_bits):
            bitlin.append((ai, red, uv))
    print(f"main-comp atom degree hist: {dict(sorted(deg_hist.items()))}")
    print(f"main-comp huge-const atoms: {len(huge)}")
    print(f"main-comp atoms linear purely in bits: {len(bitlin)}")

    def fmt(red):
        parts = []
        for m, c in sorted(red.items(), key=lambda t: (len(t[0]), t[0])):
            if m == (): parts.append(str(c))
            else: parts.append(f"{c}*" + "*".join(f"x{v}" for v in m))
        return " + ".join(parts)

    print("\n=== sample huge-const atoms in main component ===")
    for ai, red, uv in huge[:8]:
        s = fmt(red)
        print("  ", s[:200])
    print("\n=== sample bit-linear atoms (if any) ===")
    for ai, red, uv in bitlin[:8]:
        print("  ", fmt(red)[:200])

    # how many bits does each huge atom touch, and their gating bit
    json.dump({
        "main_vars": sorted(main_vars),
        "main_bits": main_bits,
        "n_huge": len(huge),
        "n_bitlin": len(bitlin),
    }, open('solve_lab/main_comp.json', 'w'))
    print("\nwrote main_comp.json")

if __name__ == '__main__':
    main()
