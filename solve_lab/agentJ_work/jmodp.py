#!/usr/bin/env python3
"""Mod-p reduced model.

Forward-propagate the definer DAG entirely in GF(p) (handles contribute p*h == 0,
so they drop out) and evaluate every constraint atom mod p.  The reduced system
is {constraint_i == 0 (mod p)} in the 13 parameters + 2 booleans.
"""
import os, pickle, sys
import jengine as E
import jman as J

P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
NV = E.NV
polys = E.polys
definer, order, FREE = J.definer, J.order, J.FREE
CONS = sorted(set(range(E.NA)) - set(definer.values()))

# precompute per-var (inv coef mod p, rest monomials)
EVP = {}
for v, i in definer.items():
    p = polys[i]
    c = p[(v,)]
    rest = tuple((k, cc % P) for k, cc in p.items() if k != (v,))
    EVP[v] = (pow(c % P, P - 2, P), rest)


def fwd_modp(val):
    """val: list of ints (residues); modifies in place, returns it."""
    for v in order:
        e = EVP.get(v)
        if e is None or v in FREE:
            continue
        ic, rest = e
        s = 0
        for k, cc in rest:
            t = cc
            for j in k:
                t = t * val[j] % P
            s += t
        val[v] = (-s) % P * ic % P
    return val


def atom_modp(i, val):
    s = 0
    for k, c in polys[i].items():
        t = c % P
        for j in k:
            t = t * val[j] % P
        s += t
    return s % P


def residues(val):
    return {i: atom_modp(i, val) for i in CONS}


def base_state(params=None):
    val = [x % P for x in J.BASE]
    if params:
        for k, v in params.items():
            val[k] = v % P
    fwd_modp(val)
    return val


if __name__ == '__main__':
    print("constraint atoms:", len(CONS))
    val = base_state()
    r = residues(val)
    bad = [i for i, x in r.items() if x]
    print("constraints nonzero mod p at base:", len(bad), bad)
    # sanity: over Z the violated set was [8583, 30271, 35890, 35892]
    PARAMS = [6418, 8778, 12553, 14623, 14853, 16742, 22152, 22162, 22649,
              24548, 30213, 31339, 33462]
    BOOLS = [2081, 24601]
    print("\nsensitivity of constraints to each parameter (+1 probe):")
    for z in PARAMS + BOOLS + [8731, 9118]:
        v2 = [x % P for x in J.BASE]
        v2[z] = (v2[z] + 1) % P
        fwd_modp(v2)
        r2 = residues(v2)
        moved = [i for i in CONS if r2[i] != r[i]]
        print(f"  x_{z}: moves {len(moved)} constraints  {moved[:14]}")
