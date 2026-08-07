#!/usr/bin/env python3
"""Handle repair: after propagation, absorb residuals into free variables that
appear linearly in a violated atom (and are not used elsewhere in a way that
matters -- we simply re-propagate and re-score, keeping only improvements)."""
import os, pickle, sys
import jengine as E
import jman as J

polys = E.polys
NV = E.NV


def linear_free_slots(i, free):
    """vars in atom i that appear linearly with some coef and not in a quadratic
    monomial."""
    p = polys[i]
    higher = set()
    for k in p:
        if len(k) >= 2:
            higher.update(k)
    out = []
    for k, c in p.items():
        if len(k) == 1 and k[0] not in higher and k[0] in free:
            out.append((k[0], c))
    return out


def atomval(i, val):
    s = 0
    for k, c in polys[i].items():
        t = c
        for j in k:
            t *= val[j]
        s += t
    return s


def repair(val, rounds=8, verbose=True):
    val = list(val)
    for r in range(rounds):
        E.forward(val, J.order, J.ev, J.definer, J.FREE)
        s, fails, av = E.score(val)
        nz = [i for i, x in enumerate(av) if x]
        if verbose:
            print(f"   repair r{r}: score={s} nz={len(nz)}")
        if not nz:
            return val, s, fails
        changed = False
        for i in nz:
            v = atomval(i, val)
            if v == 0:
                continue
            for (z, c) in linear_free_slots(i, J.FREE):
                if v % c == 0:
                    val[z] -= v // c
                    changed = True
                    break
        if not changed:
            break
    E.forward(val, J.order, J.ev, J.definer, J.FREE)
    s, fails, av = E.score(val)
    return val, s, fails


if __name__ == '__main__':
    val, s, fails = repair(J.BASE)
    print("after repair:", s, fails)
