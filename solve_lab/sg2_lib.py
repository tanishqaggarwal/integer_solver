#!/usr/bin/env python3
"""sg2 exploration library: load gates + atoms, provide lookups."""
import json, re
from collections import defaultdict
p = 2**256 - 2**32 - 977
NVARS = 38748
HERE = __file__.rsplit('/', 1)[0]

# ---- gates ----
gates = []
with open(HERE + '/atoms/gates.jsonl') as f:
    for line in f:
        d = json.loads(line)
        gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_out = set(t for t, _, _ in gates)
freeinp = set(v for v in range(NVARS) if v not in gate_out)

# gate defs indexed by target
gates_by_target = defaultdict(list)
for gi, (t, rhs, vids) in enumerate(gates):
    gates_by_target[t].append(gi)

# ---- atoms ----
def load_atoms_full():
    """Return list of dicts with keys poly(dict), repr, n_eq, eqs, idx."""
    out = []
    with open(HERE + '/atoms/poly_atoms.jsonl') as f:
        for i, line in enumerate(f):
            d = json.loads(line)
            poly = {tuple(m): c for m, c in d['poly']}
            out.append({'idx': i, 'poly': poly, 'repr': d.get('repr', ''),
                        'n_eq': d.get('n_eq', 0), 'eqs': d.get('eqs', [])})
    return out

def atom_vars(poly):
    s = set()
    for m in poly:
        s.update(m)
    return s

def gdef(t):
    """Print gate definition(s) for target t."""
    for gi in gates_by_target.get(t, []):
        _, rhs, vids = gates[gi]
        print(f"  gate {gi}: x_{t} = {rhs}   vids={list(vids)}")
    if t in freeinp:
        print(f"  x_{t} is a FREE input (no gate)")

if __name__ == '__main__':
    print(f"gates={len(gates)} targets={len(gate_out)} free={len(freeinp)}")
