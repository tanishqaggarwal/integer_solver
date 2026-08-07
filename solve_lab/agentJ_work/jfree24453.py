#!/usr/bin/env python3
"""a39032 = (x_24453) - C sits in exactly ONE equation (eq27494).

So breaking that single pin costs weight 1 (score 39032 if nothing else breaks) and
buys x_24453 as a free knob.  Does that knob reach the residual constraints?
"""
import sys, itertools
from collections import defaultdict
import jengine as E, jman as J, jmodp as MP, jsolve2 as S

P = MP.P
A_PIN = 39032
V = 24453

definer = dict(J.definer)
rm = [v for v, i in definer.items() if i == A_PIN]
print("pin atom a%d defines: %s" % (A_PIN, rm))
for v in rm:
    del definer[v]
order, cyc = E.topo(definer)
assert not cyc
FREE = set(range(E.NV)) - set(definer)
EVP = {}
for v, i in definer.items():
    p = E.polys[i]
    EVP[v] = (pow(p[(v,)] % P, P - 2, P), tuple((k, cc % P) for k, cc in p.items() if k != (v,)))
CONS = sorted(set(range(len(E.polys))) - set(definer.values()) - {A_PIN})
print("constraint atoms (excluding the broken pin):", len(CONS))


def fwd(val):
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


def residues(val):
    return {i: MP.atom_modp(i, val) for i in CONS}


for b1, b2 in [(1, 1), (1, 0), (0, 1), (0, 0)]:
    base, _ = S.branch(b1, b2)
    val = list(base)
    fwd(val)
    r0 = residues(val)
    bad0 = sorted(i for i, x in r0.items() if x)
    # does x_24453 move anything?
    v2 = list(val); v2[V] = (v2[V] + 1) % P
    fwd(v2)
    r1 = residues(v2)
    moved = [i for i in CONS if r1[i] != r0[i]]
    print(f"\nbranch ({b1},{b2}): violated {bad0}")
    print(f"   x_24453 moves {len(moved)} constraints: {moved[:12]}")
    if not moved:
        continue
    # try to solve each violated constraint with this knob (affine probe)
    val2 = list(val)
    for c in bad0:
        old = val2[V]
        val2[V] = 0; fwd(val2); a0 = MP.atom_modp(c, val2)
        val2[V] = 1; fwd(val2); a1 = MP.atom_modp(c, val2)
        g = (a1 - a0) % P
        if g:
            val2[V] = (-a0) * pow(g, P - 2, P) % P
            fwd(val2)
            r = residues(val2)
            nb = sorted(i for i, x in r.items() if x)
            print(f"   solving a{c} with x_24453 -> violated {len(nb)}: {nb[:10]}")
        else:
            val2[V] = old; fwd(val2)
            print(f"   a{c}: zero gradient in x_24453")
