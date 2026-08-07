#!/usr/bin/env python3
"""agent AF, step 3: definition DAG, constants, aliases, and the modulus P."""
import sys, os, pickle
from collections import Counter, defaultdict
sys.setrecursionlimit(100000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af2_atoms import shape_of, varsof

def main():
    D = pickle.load(open(os.path.join(HERE, 'af_atoms.pkl'), 'rb'))
    atoms = D['atoms']; eq_atoms = D['eq_atoms']

    # which equations does each atom appear in
    a2e = defaultdict(set)
    for i, la in enumerate(eq_atoms):
        for c, a in la:
            a2e[a].add(i)

    defs = defaultdict(list)      # var -> list of (atom_id, rhsAST)
    consts = {}                   # var -> int   (from (V - K) atoms)
    for aid, a in enumerate(atoms):
        if a[0] == '-' and a[1][0] == 'v':
            v = a[1][1]
            defs[v].append((aid, a[2]))
            if a[2][0] == 'c':
                consts.setdefault(v, set()).add(a[2][1])
    print('vars with at least one (V - RHS) atom: %d' % len(defs))
    print('vars pinned to a literal constant: %d' % len(consts))
    bad = {v: s for v, s in consts.items() if len(s) > 1}
    print('  vars with conflicting constant pins: %d' % len(bad))
    consts = {v: next(iter(s)) for v, s in consts.items()}

    big = sorted(((val, v) for v, val in consts.items() if val > 2**200), reverse=True)
    print('constants > 2^200: %d' % len(big))
    cc = Counter(val for val, v in big)
    for val, k in cc.most_common(6):
        print('  x%d  ... %d vars share value  bits=%d  val=%d' % (
            [v for vv, v in big if vv == val][0], k, val.bit_length(), val))
    pickle.dump({'defs': dict(defs), 'consts': consts, 'a2e': dict(a2e)},
                open(os.path.join(HERE, 'af_defs.pkl'), 'wb'), 2)

if __name__ == '__main__':
    main()
