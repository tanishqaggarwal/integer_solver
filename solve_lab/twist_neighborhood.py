#!/usr/bin/env python3
"""Extract the twist-neighborhood subsystem ORIENTATION-FREE. For x_18274 and
x_17728 (the 233-side twist targets) and x_9770,x_3183 (the 22-side), dump every
atom containing them in full polynomial form. This shows the TRUE defining
constraints (residue loads / products), independent of any forward-eval choice.
Then BFS the constraint graph from the twist vars to see how large the truly
coupled subsystem is."""
import json
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS

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

    KEY = [18274, 17728, 9770, 3183, 8821]
    for t in KEY:
        ats = var_atoms[t]
        print(f"\n===== x_{t}: appears in {len(ats)} atoms; best={val(t)} =====")
        for a in ats:
            poly = atoms[a]
            s = 0
            for m, c in poly.items():
                tt = c
                for x in m: tt *= val(x)
                s += tt
            nv = len(atom_vars(poly))
            print(f"  a{a} [nv={nv} resid@best={s}]: {fmt(poly, control)[:200]}")

    # BFS coupling from twist vars, but STOP at 'huge' residue atoms boundary to
    # measure the tightly-coupled core size.
    seen = set(KEY); q = deque(KEY); depth = {v:0 for v in KEY}
    coupled_atoms = set()
    while q:
        v = q.popleft()
        if depth[v] >= 3: continue
        for a in var_atoms[v]:
            coupled_atoms.add(a)
            for u in atom_vars(atoms[a]):
                if u not in seen:
                    seen.add(u); depth[u] = depth[v]+1; q.append(u)
    print(f"\n\nBFS depth<=3 from twist vars: {len(seen)} vars, {len(coupled_atoms)} atoms")

if __name__ == '__main__':
    main()
