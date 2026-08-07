#!/usr/bin/env python3
"""Agent Z, TASK 3 (decisive): does ANY equation constrain |S| ?

Method (named: 'booleanity-reduced affine elimination'):
  1. Find every variable that carries a booleanity atom  c*(x - x^2)  anywhere.
  2. Expand every equation's linear form L into a polynomial, then reduce it
     modulo  x^2 -> x  for every boolean variable.
  3. Keep the equations whose reduced polynomial has total degree <= 1.
     Each is an exact linear constraint valid on the boolean locus.
  4. Gaussian-eliminate the NON-selector variables out of that linear system
     (over a large prime field, then re-verified over Q on the survivors).
     Anything that survives is a linear constraint on the selector vector --
     i.e. exactly the place a cardinality constraint would live.
"""
import os, sys, json, pickle, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zparse import parse, varset, reduce_L, atoms_of
from zatoms import poly, pkey

HERE = os.path.dirname(os.path.abspath(__file__))
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')
sel = set(json.load(open(os.path.join(HERE, 'zsel.json')))['selectors'])
PRIME = (1 << 61) - 1

def main():
    sys.setrecursionlimit(100000)
    lines = [l.strip().rsplit('=', 1)[0] for l in open(EQ) if l.strip()]

    # ---- pass 1: expand every equation, record polynomial, find boolean vars
    print("pass 1: expanding all equations")
    eqpoly = []
    boolvars = set()
    for i, lhs in enumerate(lines):
        E = parse(lhs)
        L, _ = reduce_L(E)
        p = {}
        for c, a in atoms_of(L):
            pa = poly(a)
            for m, cc in pa.items():
                p[m] = p.get(m, 0) + c * cc
            # booleanity detection on the ATOM (c*(x - x^2))
            if len(pa) == 2:
                items = sorted(pa.items(), key=lambda kv: len(kv[0]))
                (m1, c1), (m2, c2) = items
                if len(m1) == 1 and len(m2) == 2 and m2 == (m1[0], m1[0]) and c1 == -c2:
                    boolvars.add(m1[0])
        p = {m: c for m, c in p.items() if c}
        eqpoly.append(p)
        if i % 5000 == 0:
            print("  ...", i, flush=True)
    print("variables carrying a booleanity atom:", len(boolvars))
    print("  of which selectors:", len(boolvars & sel), "/", len(sel))
    print("  non-selector boolean vars:", len(boolvars - sel))

    # ---- pass 2: reduce mod x^2 -> x for boolean vars, keep degree<=1 equations
    def red(mon):
        out = []
        cnt = collections.Counter(mon)
        for v, k in cnt.items():
            if v in boolvars:
                out.append(v)          # x^k -> x
            else:
                out.extend([v] * k)
        return tuple(sorted(out))

    lin_rows = []
    degcount = collections.Counter()
    for i, p in enumerate(eqpoly):
        q = {}
        for m, c in p.items():
            mm = red(m)
            q[mm] = q.get(mm, 0) + c
        q = {m: c for m, c in q.items() if c}
        d = max((len(m) for m in q), default=0)
        degcount[d] += 1
        if d <= 1:
            lin_rows.append((i, q))
    print("degree histogram after booleanity reduction:", sorted(degcount.items()))
    print("linear (deg<=1) equations after reduction:", len(lin_rows))

    pickle.dump({'boolvars': boolvars, 'lin_rows': lin_rows},
                open(os.path.join(HERE, 'zlin.pkl'), 'wb'))
    print("saved zlin.pkl")
    return

    # ---- pass 3: eliminate non-selector variables
    allvars = set()
    for i, q in lin_rows:
        for m in q:
            if m:
                allvars.add(m[0])
    nonsel = sorted(allvars - sel)
    selin = sorted(allvars & sel)
    print("vars in the linear system:", len(allvars), " non-selector:", len(nonsel),
          " selectors:", len(selin))

    # column order: non-selectors first (to be eliminated), then selectors, then const
    order = {v: j for j, v in enumerate(nonsel)}
    NS = len(nonsel)
    for j, v in enumerate(selin):
        order[v] = NS + j
    CONST = NS + len(selin)

    rows = []
    for i, q in lin_rows:
        r = {}
        for m, c in q.items():
            j = CONST if not m else order[m[0]]
            r[j] = (r.get(j, 0) + c) % PRIME
        r = {j: c for j, c in r.items() if c}
        if r:
            rows.append((i, r))
    print("rows:", len(rows))

    pivots = {}          # col -> (rowdict, origin eq)
    survivors = []
    for cnt, (i, r) in enumerate(rows):
        r = dict(r)
        prov = [i]
        while True:
            cols = [j for j in r if j < NS]
            if not cols:
                break
            j = min(cols)
            if j not in pivots:
                inv = pow(r[j], PRIME - 2, PRIME)
                r = {k: (v * inv) % PRIME for k, v in r.items()}
                pivots[j] = (r, prov)
                r = None
                break
            pr, pprov = pivots[j]
            f = r[j]
            for k, v in pr.items():
                nv = (r.get(k, 0) - f * v) % PRIME
                if nv:
                    r[k] = nv
                elif k in r:
                    del r[k]
            prov = prov + pprov
        if r is not None and r:
            survivors.append((sorted(r.items()), prov))
        if cnt % 500 == 0:
            print("   elim", cnt, "pivots", len(pivots), "survivors", len(survivors), flush=True)

    print()
    print("=" * 78)
    print("PIVOTS on non-selector columns:", len(pivots))
    print("ROWS SURVIVING elimination of every non-selector variable:", len(survivors))
    nontriv = []
    for r, prov in survivors:
        cols = [j for j, v in r if j != CONST]
        if cols:
            nontriv.append((r, prov))
    print("  ... of which mention at least one selector:", len(nontriv))
    for r, prov in survivors[:40]:
        desc = []
        for j, v in r:
            vv = v if v < PRIME // 2 else v - PRIME
            if j == CONST:
                desc.append("%+d" % vv)
            else:
                nm = ("s%d" % selin[j - NS]) if j >= NS else ("x%d" % nonsel[j])
                desc.append("%+d*%s" % (vv, nm))
        print("   SURVIVOR from eqs", sorted(set(prov))[:6], ":", " ".join(desc)[:300])

    pickle.dump({'boolvars': boolvars, 'lin_rows': lin_rows,
                 'survivors': survivors, 'nonsel': nonsel, 'selin': selin},
                open(os.path.join(HERE, 'zaffine.pkl'), 'wb'))
    print("saved zaffine.pkl")

if __name__ == '__main__':
    main()
