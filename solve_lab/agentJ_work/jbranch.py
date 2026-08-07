#!/usr/bin/env python3
"""For each boolean branch (x_2081, x_24601, x_4287, ...) run a generic
Gauss-Seidel over the 13 parameters and report the remaining violated
constraints mod p."""
import os, sys, itertools, random
import jengine as E, jman as J, jmodp as MP

P = MP.P
PARAMS = [6418, 8778, 12553, 14623, 14853, 16742, 22152, 22162, 22649,
          24548, 30213, 31339, 33462]


def probe_lin(val, c, z):
    """(a,b) such that atom c == a*z + b mod p (assumes affine)."""
    old = val[z]
    val[z] = 0; MP.fwd_modp(val); r0 = MP.atom_modp(c, val)
    val[z] = 1; MP.fwd_modp(val); r1 = MP.atom_modp(c, val)
    val[z] = old; MP.fwd_modp(val)
    return (r1 - r0) % P, r0


def sweep(val, rounds=15, verbose=False):
    for it in range(rounds):
        r = MP.residues(val)
        bad = [i for i, x in r.items() if x]
        if not bad:
            return bad
        progress = False
        for c in bad:
            for z in PARAMS:
                a, b = probe_lin(val, c, z)
                if a:
                    val[z] = (-b) * pow(a, P - 2, P) % P
                    MP.fwd_modp(val)
                    progress = True
                    break
        r = MP.residues(val)
        nb = [i for i, x in r.items() if x]
        if verbose:
            print(f"    round {it}: {len(nb)} violated {nb[:12]}")
        if set(nb) == set(bad) and not progress:
            break
        bad = nb
    return bad


if __name__ == '__main__':
    combos = list(itertools.product([0, 1], repeat=3))
    for b1, b2, b3 in combos:
        val = [x % P for x in J.BASE]
        val[2081] = b1; val[24601] = b2; val[4287] = b3
        MP.fwd_modp(val)
        r = MP.residues(val)
        start = sorted(i for i, x in r.items() if x)
        bad = sweep(val)
        print(f"x2081={b1} x24601={b2} x4287={b3}: start {len(start)} -> after sweep {len(bad)}: {sorted(bad)}")
