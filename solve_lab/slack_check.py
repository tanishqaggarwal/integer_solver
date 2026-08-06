#!/usr/bin/env python3
"""Are x_26977 and x_9982 slack vars (appear in only one atom each)? If so we can
set them to satisfy atoms 1817 / 30378 directly. Check incidence across all atoms;
also which of the twist wires appear where."""
import json
from propagate import load_atoms, atom_vars

def main():
    A = load_atoms()
    watch = [26977, 9982, 9770, 18274, 3183, 17728]
    inc = {w: [] for w in watch}
    for a, poly in enumerate(A):
        va = atom_vars(poly)
        for w in watch:
            if w in va: inc[w].append(a)
    for w in watch:
        print(f"x_{w}: appears in {len(inc[w])} atoms: {inc[w][:12]}")
    # detail: for the slack candidates, show the atoms
    for w in (26977, 9982):
        print(f"\n--- atoms containing x_{w} ---")
        for a in inc[w]:
            terms = []
            for m, c in A[a].items():
                terms.append(f"{c:+d}*{'*'.join('x_'+str(x) for x in m) if m else '1'}")
            print(f"  atom {a}: {' '.join(terms)[:200]}")

if __name__ == '__main__':
    main()
