#!/usr/bin/env python3
"""Reverse-mode AD over GF(q): sensitivity of the residual constraints to EVERY
variable, including defined ones (i.e. what you would gain by breaking that
variable's pin).  One pass per target instead of one per knob.

Breaking the definer atom of a variable v costs |eqs(definer(v))| equations (those
become the nonzero-atom support) and buys v as a free knob.  We want a cheap set of
pins whose gradients span the residual.
"""
import sys, pickle, os
from collections import defaultdict
import jengine as E, jman as J, jmodp as MP, jsolve2 as S
import jdist as DI

P = MP.P
definer, order = J.definer, J.order
polys = E.polys
pos = {v: k for k, v in enumerate(order)}


def datom_dvar(i, v, val):
    """d(atom i)/d(x_v) mod q"""
    g = 0
    for k, c in polys[i].items():
        if v not in k:
            continue
        if len(k) == 1:
            g += c
        else:
            a, b = k
            if a == b:
                g += 2 * c * val[a]
            else:
                g += c * val[b if a == v else a]
    return g % P


def reverse_grad(target, val):
    """g[v] = d(atom target)/d(x_v) for every variable v, through the definer DAG."""
    g = defaultdict(int)
    for v in E.varsof[target]:
        g[v] = (g[v] + datom_dvar(target, v, val)) % P
    for v in reversed(order):
        gv = g.get(v, 0)
        if not gv:
            continue
        i = definer.get(v)
        if i is None:
            continue
        c = polys[i][(v,)] % P
        inv = pow(c, P - 2, P)
        for u in E.varsof[i]:
            if u == v:
                continue
            d = datom_dvar(i, u, val)
            if d:
                g[u] = (g[u] - gv * d % P * inv) % P
    return {k: x for k, x in g.items() if x}


if __name__ == '__main__':
    b1, b2 = 1, 1
    val, bad = S.branch(b1, b2)
    print(f"branch ({b1},{b2}) violated {bad}")
    grads = {}
    for c in bad:
        grads[c] = reverse_grad(c, val)
        print(f"  a{c}: {len(grads[c])} variables have nonzero sensitivity")

    movers = set()
    for c in bad:
        movers |= set(grads[c])
    print("variables moving at least one residual constraint:", len(movers))

    # cost of breaking each: |eqs(definer(v))|
    rows = []
    for v in movers:
        i = definer.get(v)
        if i is None:
            rows.append((0, v, None))          # already free
        else:
            rows.append((len(DI.A2E[i]), v, i))
    rows.sort()
    print("\ncheapest pins to break (cost = #equations of the definer atom):")
    for cost, v, i in rows[:30]:
        vec = tuple(grads[c].get(v, 0) != 0 for c in bad)
        print(f"   cost {cost:3d}  x_{v:<6} definer a{i}  moves {vec}")
    pickle.dump({'grads': grads, 'bad': bad, 'rows': rows},
                open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'jrev.pkl'), 'wb'))
