#!/usr/bin/env python3
"""After propagation, characterize the residual 'hard core':
 - unassigned variables (which are booleans)
 - unsolved atoms (>=2 unknowns after substitution)
 - connected components of the residual constraint graph
"""
import json, time
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, substitute

NVARS = 38748

def main():
    atoms = load_atoms()
    bool_vars = set()
    for poly in atoms:
        d = {tuple(m): c for m, c in poly.items()}
        vs = atom_vars(poly)
        if len(vs) == 1 and len(d) == 2:
            v = next(iter(vs))
            if d.get((v,)) == 1 and d.get((v, v)) == -1:
                bool_vars.add(v)
            elif d.get((v,)) == -1 and d.get((v, v)) == 1:
                bool_vars.add(v)
    print(f"boolean vars: {len(bool_vars)}")

    part = json.load(open('solve_lab/partial_assignment.json'))
    val = [None] * NVARS
    for k, x in part.items():
        val[int(k[2:])] = x
    n_assigned = sum(1 for x in val if x is not None)
    print(f"assigned from propagation: {n_assigned}")

    # unsolved atoms
    unsolved = []
    for ai, poly in enumerate(atoms):
        red = substitute(poly, val)
        uv = atom_vars(red)
        if len(uv) >= 2:
            unsolved.append((ai, uv))
    print(f"unsolved atoms (>=2 unknowns): {len(unsolved)}")

    # union-find over unassigned vars connected by unsolved atoms
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb

    for ai, uv in unsolved:
        uv = list(uv)
        for v in uv:
            find(v)
        for i in range(1, len(uv)):
            union(uv[0], uv[i])

    comp = defaultdict(set)
    for v in list(parent.keys()):
        comp[find(v)].add(v)
    sizes = sorted((len(s) for s in comp.values()), reverse=True)
    print(f"residual components: {len(comp)}")
    print(f"component size histogram (top 20): {sizes[:20]}")
    from collections import Counter
    sc = Counter(sizes)
    print(f"size distribution: {dict(sorted(sc.items()))}")

    # free bits (unassigned booleans)
    free_bits = [v for v in bool_vars if val[v] is None]
    print(f"free (unassigned) boolean vars: {len(free_bits)}")
    # unassigned non-bit vars
    unassigned = [v for v in range(NVARS) if val[v] is None]
    print(f"total unassigned: {len(unassigned)}")
    print(f"unassigned that are NOT booleans: {len(unassigned)-len(free_bits)}")

    # per-component bit count
    comp_bits = []
    for root, s in comp.items():
        nb = sum(1 for v in s if v in bool_vars)
        comp_bits.append((len(s), nb))
    comp_bits.sort(reverse=True)
    print("largest components (size, #bits):", comp_bits[:15])

    json.dump({
        "n_bool": len(bool_vars),
        "n_assigned": n_assigned,
        "n_unsolved_atoms": len(unsolved),
        "n_components": len(comp),
        "comp_size_hist": dict(sorted(sc.items())),
        "n_free_bits": len(free_bits),
        "largest_comps": comp_bits[:30],
    }, open('solve_lab/core_analysis.json', 'w'), indent=1)
    print("wrote core_analysis.json")

if __name__ == '__main__':
    main()
