#!/usr/bin/env python3
"""Solve the linear chain of the reduced mod-p system by Gauss-Seidel, then
report the three EC constraints."""
import os, sys, random
import jengine as E, jman as J, jmodp as MP

P = MP.P
EC = [20407, 20409, 31575]
# (constraint, knob) pairs: knob appears with degree 1 in the constraint
PAIRS = [(3895, 6418), (3897, 12553), (32257, 22152), (32259, 33462),
         (30271, 14853), (8583, 24548), (22688, 14623), (26603, 31339),
         (34370, 8778), (27640, 16742), (2694, 22649),
         (31571, 22162), (731, 30213)]
ALLC = [c for c, _ in PAIRS] + EC


def solve_pair(val, c, z):
    """set val[z] so that atom c vanishes mod p (assumes affine in z)."""
    old = val[z]
    val[z] = 0; MP.fwd_modp(val); r0 = MP.atom_modp(c, val)
    val[z] = 1; MP.fwd_modp(val); r1 = MP.atom_modp(c, val)
    a = (r1 - r0) % P
    if a == 0:
        val[z] = old; MP.fwd_modp(val)
        return None
    znew = (-r0) * pow(a, P - 2, P) % P
    val[z] = znew; MP.fwd_modp(val)
    return znew


def report(val, tag):
    r = {c: MP.atom_modp(c, val) for c in ALLC}
    bad = [c for c in ALLC if r[c]]
    print(f"  {tag}: violated {bad}")
    return bad


if __name__ == '__main__':
    val = [x % P for x in J.BASE]
    MP.fwd_modp(val)
    report(val, 'start')
    for it in range(12):
        for c, z in PAIRS:
            solve_pair(val, c, z)
        bad = report(val, f'sweep {it}')
        if not bad:
            break
    print("\nfull constraint check mod p:")
    r = MP.residues(val)
    print("violated:", sorted(i for i, x in r.items() if x))
    import pickle
    pickle.dump(val, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'jchain_val.pkl'), 'wb'))
