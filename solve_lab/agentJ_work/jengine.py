#!/usr/bin/env python3
"""Generative engine: definer map -> topological order -> forward evaluation.

definer(v) = the atom whose text starts with '(x_v)-' and for which v is a legal
def target.  Ties broken by first occurrence; the losers become constraints.
"""
import os, pickle, json, sys
from collections import defaultdict, deque, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
M = pickle.load(open(os.path.join(HERE, 'jmodel2.pkl'), 'rb'))
P = pickle.load(open(os.path.join(HERE, 'jpoly.pkl'), 'rb'))
LEAD = pickle.load(open(os.path.join(HERE, 'jlead.pkl'), 'rb'))
eqs, atoms, polys, defcands = M['eqs'], M['atoms'], P['polys'], P['defcands']
NA = len(polys); NV = 38748

varsof = []
for p in polys:
    s = set()
    for k in p:
        s.update(k)
    varsof.append(s)


def build_definer(lead=LEAD):
    definer = {}
    for i, v in enumerate(lead):
        if v is None:
            continue
        if v not in definer:
            definer[v] = i
    return definer


def topo(definer):
    """Return (order, cyclic) — order is list of vars in evaluation order."""
    indeg = Counter()
    succ = defaultdict(list)
    for v, i in definer.items():
        for w in varsof[i]:
            if w != v:
                succ[w].append(v)
                indeg[v] += 1
    q = deque([v for v in range(NV) if indeg[v] == 0])
    order = []
    while q:
        v = q.popleft()
        order.append(v)
        for w in succ[v]:
            indeg[w] -= 1
            if indeg[w] == 0:
                q.append(w)
    cyc = [v for v in range(NV) if indeg[v] > 0]
    return order, cyc


def make_eval(definer):
    """Precompute per-var (coef_of_v, rest_poly) so v = -rest/coef."""
    ev = {}
    for v, i in definer.items():
        p = polys[i]
        c = p[(v,)]
        rest = {k: cc for k, cc in p.items() if k != (v,)}
        ev[v] = (c, tuple(rest.items()))
    return ev


def forward(val, order, ev, definer, free):
    """Evaluate all defined vars in topological order.  Integer division must be
    exact; returns list of vars whose division was inexact."""
    bad = []
    for v in order:
        if v in free:
            continue
        e = ev.get(v)
        if e is None:
            continue
        c, rest = e
        s = 0
        for k, cc in rest:
            t = cc
            for j in k:
                t *= val[j]
            s += t
        q, r = divmod(-s, c)
        if r != 0:
            bad.append(v)
            q = -s // c
        val[v] = q
    return bad


def atomvals(val):
    out = []
    for p in polys:
        s = 0
        for k, c in p.items():
            t = c
            for j in k:
                t *= val[j]
            s += t
        out.append(s)
    return out


def score(val):
    av = atomvals(val)
    fails = []
    for e in eqs:
        s = 0
        for c, j in e['terms']:
            s += c * av[j]
        if s != 0:
            fails.append(e['i'])
    return len(eqs) - len(fails), fails, av


def load(path):
    d = json.load(open(path))
    val = [0] * NV
    for k, v in d.items():
        val[int(k[2:]) if k.startswith('x_') else int(k)] = int(v)
    return val


def save(val, path):
    json.dump({f"x_{i}": val[i] for i in range(NV) if val[i] != 0}, open(path, 'w'))


if __name__ == '__main__':
    definer = build_definer()
    order, cyc = topo(definer)
    print("definer covers", len(definer), "vars; topo order", len(order), "cyclic", len(cyc))
    ev = make_eval(definer)
    val = load(sys.argv[1] if len(sys.argv) > 1 else
               os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json'))
    ref = list(val)
    free = set(range(NV)) - set(definer)
    print("free vars:", len(free))
    bad = forward(val, order, ev, definer, free)
    print("inexact divisions:", len(bad))
    diff = [i for i in range(NV) if val[i] != ref[i]]
    print("vars changed by re-derivation:", len(diff), diff[:20])
    s, fails, av = score(val)
    print("score after re-derivation:", s, fails[:12])
