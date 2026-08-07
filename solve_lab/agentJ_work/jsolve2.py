#!/usr/bin/env python3
"""Fast targeted chain sweep per boolean branch (mod p)."""
import sys, itertools
import jengine as E, jman as J, jmodp as MP
P = MP.P
PAIRS = [(3895, 6418), (3899, 6418), (3897, 12553), (3901, 12553),
         (32257, 22152), (1162, 22152), (32259, 33462), (1164, 33462),
         (30271, 14853), (8583, 24548), (22688, 14623), (26603, 31339),
         (34370, 8778), (27640, 16742), (2694, 22649),
         (31571, 22162), (31567, 22162), (731, 30213), (31569, 30213)]


def solve_pair(val, c, z):
    old = val[z]
    val[z] = 0; MP.fwd_modp(val); r0 = MP.atom_modp(c, val)
    val[z] = 1; MP.fwd_modp(val); r1 = MP.atom_modp(c, val)
    a = (r1 - r0) % P
    if a == 0:
        val[z] = old; MP.fwd_modp(val); return False
    val[z] = (-r0) * pow(a, P - 2, P) % P
    MP.fwd_modp(val)
    return True


def branch(b1, b2, b3=0, extra=None, rounds=8, verbose=False):
    val = [x % P for x in J.BASE]
    val[2081] = b1; val[24601] = b2; val[4287] = b3
    val[9118] = 0; val[8731] = 0
    if extra:
        for k, v in extra.items():
            val[k] = v % P
    MP.fwd_modp(val)
    bad = None
    for it in range(rounds):
        r = MP.residues(val)
        bad = set(i for i, x in r.items() if x)
        if not bad:
            break
        for c, z in PAIRS:
            if MP.atom_modp(c, val) != 0:
                solve_pair(val, c, z)
        r = MP.residues(val)
        nb = set(i for i, x in r.items() if x)
        if verbose:
            print('   round', it, len(nb), sorted(nb)[:14])
        if nb == bad:
            bad = nb
            break
        bad = nb
    r = MP.residues(val)
    return val, sorted(i for i, x in r.items() if x)


if __name__ == '__main__':
    for b1, b2 in itertools.product([0, 1], repeat=2):
        val, bad = branch(b1, b2)
        print(f"b1={b1} b2={b2}: remaining violated mod p = {len(bad)}  {bad}")
