#!/usr/bin/env python3
"""Find CONFINED cheap pins.

Breaking the definer atom of v costs |eqs(definer(v))| equations and frees v.  The knob
is only useful if its influence is CONFINED: perturbing v must disturb no constraint
outside the residual (otherwise we are back to the frozen situation).  x_24453 has this
property -- it moves exactly the 3 residual constraints and nothing else.

Search all defined variables that move the residual, cheapest pin first, and report
which constraints each one disturbs.
"""
import sys, pickle, os
from collections import defaultdict
import jengine as E, jman as J, jmodp as MP, jsolve2 as S
import jdist as DI

P = MP.P
definer, order = J.definer, J.order
polys = E.polys
HERE = os.path.dirname(os.path.abspath(__file__))

R = pickle.load(open(os.path.join(HERE, 'jrev.pkl'), 'rb'))
grads, bad = R['grads'], R['bad']
movers = set()
for c in bad:
    movers |= set(grads[c])

cands = []
for v in movers:
    i = definer.get(v)
    if i is not None:
        cands.append((len(DI.A2E[i]), v, i))
cands.sort()
print(f"defined movers (breakable pins): {len(cands)}; cheapest costs "
      f"{[c[0] for c in cands[:20]]}")

base, _ = S.branch(1, 1)
badset = set(bad)


def released_effect(v, pin):
    """release pin (v becomes free), perturb v, report disturbed constraints."""
    d2 = dict(definer)
    for w in [w for w, i in d2.items() if i == pin]:
        del d2[w]
    order2, cyc = E.topo(d2)
    if cyc:
        return None
    FREE2 = set(range(E.NV)) - set(d2)
    EVP2 = {}
    for w, i in d2.items():
        p = polys[i]
        EVP2[w] = (pow(p[(w,)] % P, P - 2, P),
                   tuple((k, cc % P) for k, cc in p.items() if k != (w,)))
    CONS2 = sorted(set(range(len(polys))) - set(d2.values()) - {pin})

    def fwd(val):
        for w in order2:
            e = EVP2.get(w)
            if e is None or w in FREE2:
                continue
            ic, rest = e
            s = 0
            for k, cc in rest:
                t = cc
                for j in k:
                    t = t * val[j] % P
                s += t
            val[w] = (-s) % P * ic % P
        return val

    v0 = fwd(list(base))
    r0 = {i: MP.atom_modp(i, v0) for i in CONS2}
    v1 = list(v0); v1[v] = (v1[v] + 1) % P
    fwd(v1)
    r1 = {i: MP.atom_modp(i, v1) for i in CONS2}
    moved = [i for i in CONS2 if r1[i] != r0[i]]
    viol = sorted(i for i, x in r0.items() if x)
    return moved, viol


print("\ncost | var | pin | #constraints disturbed | confined to residual?")
good = []
for cost, v, pin in cands[:60]:
    res = released_effect(v, pin)
    if res is None:
        continue
    moved, viol = res
    conf = set(moved) <= badset
    print(f" {cost:4d}  x_{v:<6} a{pin:<6} moved={len(moved):<4} confined={conf}"
          f" {'<== USABLE' if conf and moved else ''}")
    if conf and moved:
        good.append((cost, v, pin, tuple(sorted(moved))))
print("\nCONFINED cheap pins:", good)
pickle.dump(good, open(os.path.join(HERE, 'jpins.pkl'), 'wb'))
