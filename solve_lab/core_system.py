#!/usr/bin/env python3
"""Extract and display the CORE algebraic subsystem: the atoms among the wires that
the twist reduces to. Goal: find the setter's construction identity (e.g. a quadratic
in x_18274 with constant coeffs, or a product relation linking 22-side to 233-side)
that would let us INVERT the trapdoor."""
import json
from collections import deque, defaultdict
from propagate import load_atoms, atom_vars

def main():
    atoms = load_atoms()
    best = json.load(open('best/best_partial_39019.json'))
    bv = {int(k[2:]): v for k, v in best.items()}
    control = set(json.load(open('control_bits.json')))
    def val(v): return bv.get(v, 0)
    def cls(v):
        x = abs(val(v))
        if v in control: return 'BIT'
        return '0' if x == 0 else (f'={val(v)}' if x < 10**7 else f'H{len(str(x))}')

    seed = [9770, 3183, 18274, 17728]
    # build var->atoms index
    var_atoms = defaultdict(list)
    for a, poly in enumerate(atoms):
        for v in atom_vars(poly): var_atoms[v].append(a)

    # grow core: add vars reachable via SHORT atoms (<=4 vars) from seed, up to a cap
    core = set(seed)
    frontier = deque(seed)
    LIMIT = 60
    while frontier and len(core) < LIMIT:
        v = frontier.popleft()
        # find the shortest atom defining v (fewest vars)
        cand = sorted(var_atoms[v], key=lambda a: len(atom_vars(atoms[a])))
        for a in cand[:3]:
            va = atom_vars(atoms[a])
            if len(va) <= 5:
                for x in va:
                    if x not in core and x not in control and len(core) < LIMIT:
                        core.add(x); frontier.append(x)
    print(f"core size {len(core)}")

    # print all atoms fully within core (+ allow constants/bits)
    shown = 0
    print("\n=== atoms within the core (grouped) ===")
    for a, poly in enumerate(atoms):
        va = atom_vars(poly)
        noncore = va - core - control
        if noncore: continue
        if len(poly) > 8: continue  # skip big combos here
        # format
        tt = []
        for m, c in sorted(poly.items(), key=lambda kv: (-len(kv[0]), kv[0])):
            cc = str(c) if abs(c) < 10**8 else f'H{len(str(abs(c)))}'
            vv = '*'.join(f'x{x}' + ('#' if x in control else '') for x in m) if m else '1'
            tt.append(f'{cc}*{vv}')
        print(f"  a{a}: " + ' + '.join(tt))
        shown += 1
        if shown > 60: break
    print(f"\ncore wire classes:")
    for v in sorted(core):
        print(f"  x{v}: {cls(v)}", end='   ')
        if (sorted(core).index(v)+1) % 4 == 0: print()
    print()

if __name__ == '__main__':
    main()
