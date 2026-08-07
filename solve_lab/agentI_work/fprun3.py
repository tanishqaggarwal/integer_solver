#!/usr/bin/env python3
"""Mod-p run with effective-support propagation that COLLECTS all conflicts
instead of stopping at the first, and scores the resulting mod-p state."""
import pickle, os, collections, sys, random, time, json
from model import Model, load_assign
from fp import FpEngine, P, sqrt_p

HERE = os.path.dirname(os.path.abspath(__file__))
NV = 38748


class Runner:
    def __init__(self):
        self.M = Model()
        self.E = FpEngine(self.M)
        self.wit = load_assign(os.path.join(HERE, '..', 'best',
                                            'new_instance_partial_39026.json'))
        self.witp = [w % P for w in self.wit]

    def run(self, policy='wit', seed=1, forced=None, blacklist=None,
            verbose=False):
        M, E = self.M, self.E
        rng = random.Random(seed)
        bl = set(blacklist or [])
        forced = dict(forced or {})
        val = [None] * NV
        reason = [None] * NV
        order = []
        decisions = []
        comp = E.comp; v2a = E.var2atoms
        Q = collections.deque(range(M.na))
        inq = bytearray(b'\x01') * M.na
        conflicts = []
        dead = bytearray(M.na)     # atom already reported as violated

        def reduce_atom(a):
            red = {}
            for m, c in comp[a]:
                t = c; um = []
                for x in m:
                    if val[x] is None:
                        um.append(x)
                    else:
                        t = t * val[x] % P
                if t == 0:
                    continue
                k = tuple(sorted(um))
                red[k] = (red.get(k, 0) + t) % P
            return {k: c for k, c in red.items() if c}

        def assign(u, x, r):
            val[u] = x; reason[u] = r; order.append(u)
            for b in v2a[u]:
                if not inq[b]:
                    inq[b] = 1; Q.append(b)

        def propagate():
            branch = []
            while Q:
                a = Q.popleft(); inq[a] = 0
                red = reduce_atom(a)
                vs = set()
                for k in red:
                    vs |= set(k)
                if not vs:
                    if red and not dead[a]:
                        dead[a] = 1
                        conflicts.append(a)
                    continue
                if len(vs) > 1:
                    continue
                u = next(iter(vs))
                c0 = red.get((), 0); c1 = red.get((u,), 0); c2 = red.get((u, u), 0)
                if c2 == 0:
                    assign(u, (-c0) * pow(c1, -1, P) % P, a)
                else:
                    d = (c1 * c1 - 4 * c2 * c0) % P
                    s = sqrt_p(d)
                    if s is None:
                        if not dead[a]:
                            dead[a] = 1; conflicts.append(a)
                        # pick a value anyway to keep going: u := -c1/(2c2)
                        assign(u, (-c1) * pow(2 * c2 % P, -1, P) % P, a)
                        continue
                    inv = pow(2 * c2 % P, -1, P)
                    roots = sorted({(-c1 + s) * inv % P, (-c1 - s) * inv % P})
                    if len(roots) == 1:
                        assign(u, roots[0], a)
                    else:
                        branch.append((u, a, roots))
            return branch

        while True:
            branch = propagate()
            pend = [(u, a, r) for u, a, r in branch if val[u] is None]
            if pend:
                for u, a, roots in pend:
                    if val[u] is not None:
                        continue
                    if u in forced:
                        x = forced[u] % P
                        if x not in roots:
                            x = roots[0]
                    elif policy == 'wit':
                        w = self.witp[u]; x = w if w in roots else roots[0]
                    elif policy == 'zero':
                        x = 0 if 0 in roots else roots[0]
                    else:
                        x = rng.choice(roots)
                    assign(u, x, 'dec'); decisions.append(u)
                continue
            two = collections.Counter()
            for a in range(M.na):
                red = reduce_atom(a)
                vs = set()
                for k in red:
                    vs |= set(k)
                if len(vs) == 2:
                    for x in vs:
                        two[x] += 1
            cands = [(n, v) for v, n in two.items() if v not in bl]
            if cands:
                u = max(cands)[1]
            else:
                rest = [v for v in range(NV) if val[v] is None]
                if not rest:
                    break
                r2 = [v for v in rest if v not in bl]
                u = r2[0] if r2 else rest[0]
            if u in forced:
                x = forced[u] % P
            elif policy == 'wit':
                x = self.witp[u]
            elif policy == 'zero':
                x = 0
            else:
                x = rng.randrange(P)
            assign(u, x, 'dec'); decisions.append(u)
        return {'val': val, 'reason': reason, 'order': order,
                'decisions': decisions, 'conflicts': conflicts}


def eqscore(M, E, val):
    """number of EQUATIONS whose core is nonzero mod p."""
    filled = [0 if x is None else x for x in val]
    av = [E.eval_atom(a, filled) for a in range(M.na)]
    bad = 0
    badlist = []
    for e, ts in enumerate(M.eq_terms):
        s = 0
        for c, a in ts:
            s += c * av[a]
        if s % P:
            bad += 1
            badlist.append(e)
    return bad, badlist, av


if __name__ == '__main__':
    R = Runner()
    pol = sys.argv[1] if len(sys.argv) > 1 else 'wit'
    out = R.run(pol, verbose=True)
    val = out['val']
    known = sum(1 for x in val if x is not None)
    print(f"policy={pol} known={known}/{NV} decisions={len(out['decisions'])} "
          f"atom-conflicts={len(out['conflicts'])}")
    bad, badlist, av = eqscore(R.M, R.E, val)
    print(f"equations failing mod p: {bad}/{R.M.ne}")
    print("conflict atoms:", out['conflicts'][:40])
    for a in out['conflicts'][:40]:
        print("   a%d %s" % (a, R.M.src[a][:120]))
    pickle.dump(out, open(os.path.join(HERE, f'fprun3_{pol}.pkl'), 'wb'))
