#!/usr/bin/env python3
"""The only freedom the tangent analysis cannot see: the boolean knobs.

Each boolean knob z satisfies a constraint a*z^2 + b*z = 0, so z in {0, -b/a}.
Flip each one to its non-zero value, re-propagate EXACTLY over GF(p), and report
what happens to the violated residues and to every other constraint.
"""
import sys, os, pickle
import jengine as E, jman as J, jmodp as MP
import jdiag as D

P = MP.P

if __name__ == '__main__':
    b1, b2 = int(sys.argv[1]), int(sys.argv[2])
    obj = D.build(b1, b2)
    val, bad, knobs, data, r0 = obj['val'], obj['bad'], obj['knobs'], obj['data'], obj['r0']
    base_bad = set(bad)
    print(f"branch ({b1},{b2}) violated {bad}")

    # boolean knobs and their alternative value
    alts = {}
    for z in knobs:
        for i in MP.CONS:
            if E.varsof[i] == {z}:
                p = E.polys[i]
                a = p.get((z, z), 0) % P
                b = p.get((z,), 0) % P
                if a:
                    alts[z] = (-b) * pow(a, P - 2, P) % P
                break
    print(f"boolean knobs: {len(alts)}")

    movers, neutral, costly = [], [], []
    for z, alt in sorted(alts.items()):
        cur = val[z] % P
        new = alt if cur != alt else 0
        v2 = list(val); v2[z] = new
        MP.fwd_modp(v2)
        r = MP.residues(v2)
        nb = set(i for i, x in r.items() if x)
        moved = [i for i in bad if r[i] != r0[i]]
        broke = nb - base_bad
        healed = base_bad - nb
        rec = (z, len(broke), len(healed), bool(moved))
        if moved:
            movers.append(rec)
        if not broke and not healed:
            neutral.append(z)
        if healed:
            costly.append((z, sorted(healed), len(broke)))
    print(f"bits that MOVE a violated residue: {len(movers)}")
    print(f"bits that are completely score-neutral on the constraint set: {len(neutral)}")
    print(f"bits that HEAL any violated constraint: {len(costly)}")
    for z, h, nb in costly[:20]:
        print(f"   x_{z}: heals {h}, breaks {nb}")
    print("\nmovers (bit, #broken, #healed):")
    for rec in movers[:40]:
        print("  ", rec)
