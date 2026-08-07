#!/usr/bin/env python3
"""Agent Z, TASK 4: the boolean / liveness layer.

Extract every equation all of whose variables are boolean (carry a booleanity
atom).  That subsystem contains the liveness ORs, the selector aliases and the
off-pins.  Then, for a given selector configuration sigma, unit-propagate it and
check (i) consistency, (ii) that every boolean wire lands in {0,1},
(iii) how many leaves / gates end up live.

If the subsystem is satisfiable for every sigma tried, at every weight from 0 to
256, then the liveness layer imposes no bound on |S|.
"""
import os, sys, json, pickle, collections, random
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zparse import parse, varset, reduce_L, atoms_of
from zatoms import poly

HERE = os.path.dirname(os.path.abspath(__file__))
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')
sel = set(json.load(open(os.path.join(HERE, 'zsel.json')))['selectors'])

def build():
    cache = os.path.join(HERE, 'zbool.pkl')
    if os.path.exists(cache):
        return pickle.load(open(cache, 'rb'))
    sys.setrecursionlimit(100000)
    lines = [l.strip().rsplit('=', 1)[0] for l in open(EQ) if l.strip()]
    eqp = []
    boolvars = set()
    for i, lhs in enumerate(lines):
        E = parse(lhs)
        L, _ = reduce_L(E)
        p = {}
        for c, a in atoms_of(L):
            pa = poly(a)
            for m, cc in pa.items():
                p[m] = p.get(m, 0) + c * cc
            if len(pa) == 2:
                it = sorted(pa.items(), key=lambda kv: len(kv[0]))
                (m1, c1), (m2, c2) = it
                if len(m1) == 1 and len(m2) == 2 and m2 == (m1[0], m1[0]) and c1 == -c2:
                    boolvars.add(m1[0])
        eqp.append({m: c for m, c in p.items() if c})
        if i % 5000 == 0:
            print("  parse", i, flush=True)
    D = (eqp, boolvars)
    pickle.dump(D, open(cache, 'wb'))
    return D

def main():
    eqp, boolvars = build()
    print("boolean vars:", len(boolvars), " selectors among them:", len(boolvars & sel))
    bsub = []
    for i, p in enumerate(eqp):
        vs = set()
        for m in p:
            vs |= set(m)
        if vs and vs <= boolvars:
            bsub.append((i, p, vs))
    print("equations wholly inside the boolean layer:", len(bsub))
    bvars = set()
    for i, p, vs in bsub:
        bvars |= vs
    print("boolean vars touched by that subsystem:", len(bvars),
          " selectors:", len(bvars & sel), " non-selector:", len(bvars - sel))

    # reduce x^2 -> x for boolean vars once
    def red(p):
        q = {}
        for m, c in p.items():
            mm = tuple(sorted(set(m)))
            q[mm] = q.get(mm, 0) + c
        return {m: c for m, c in q.items() if c}
    bsub = [(i, red(p), vs) for i, p, vs in bsub]
    deg = collections.Counter(max((len(m) for m in p), default=0) for i, p, vs in bsub)
    print("degree histogram of the boolean subsystem (after x^2->x):", sorted(deg.items()))

    var_eqs = collections.defaultdict(list)
    for k, (i, p, vs) in enumerate(bsub):
        for v in vs:
            var_eqs[v].append(k)

    def run(sigma_on):
        """sigma_on: set of selectors set to 1.  Returns (ok, val, nbad)."""
        val = {}
        for s in sel:
            val[s] = 1 if s in sigma_on else 0
        queue = collections.deque(range(len(bsub)))
        inq = set(queue)
        while queue:
            k = queue.popleft(); inq.discard(k)
            i, p, vs = bsub[k]
            unk = [v for v in vs if v not in val]
            if len(unk) == 0:
                tot = 0
                for m, c in p.items():
                    t = c
                    for v in m:
                        t *= val[v]
                    tot += t
                if tot != 0:
                    return False, val, ('contradiction', i)
                continue
            if len(unk) > 1:
                continue
            u = unk[0]
            a = 0; b = 0
            for m, c in p.items():
                t = c
                hit = u in m
                for v in m:
                    if v != u:
                        t *= val[v]
                if hit:
                    a += t
                else:
                    b += t
            if a == 0:
                if b != 0:
                    return False, val, ('contradiction-nolinear', i)
                continue
            x = Fraction(-b, a)
            if x not in (0, 1):
                return False, val, ('non-boolean wire', i, u, x)
            val[u] = int(x)
            for k2 in var_eqs[u]:
                if k2 not in inq:
                    queue.append(k2); inq.add(k2)
        # final full check of every fully-assigned equation
        bad = 0
        for i, p, vs in bsub:
            if all(v in val for v in vs):
                tot = 0
                for m, c in p.items():
                    t = c
                    for v in m:
                        t *= val[v]
                    tot += t
                if tot != 0:
                    bad += 1
        return True, val, bad

    rng = random.Random(11)
    S = sorted(sel)
    tests = []
    for w in [0, 1, 2, 3, 8, 17, 32, 64, 100, 128, 129, 156, 200, 240, 254, 255, 256]:
        tests.append((w, set(rng.sample(S, w))))
    for w in [128, 128, 128]:
        tests.append((w, set(rng.sample(S, w))))
    print()
    print("%5s %8s %9s %9s %s" % ("w", "ok", "resolved", "unres", "note"))
    for w, on in tests:
        ok, val, note = run(on)
        res = len([v for v in bvars if v in val])
        print("%5d %8s %9d %9d %s" % (w, ok, res, len(bvars) - res, note))

if __name__ == '__main__':
    main()
