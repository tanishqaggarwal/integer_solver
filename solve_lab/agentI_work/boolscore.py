#!/usr/bin/env python3
"""Fast objective: given the boolean selector settings, propagate mod p with
effective-support and count violated atoms / equations.  No free-input
designation (that is the slow part) -- unknowns stay unknown."""
import pickle, os, collections, sys, random, time, json
from model import Model, load_assign
from fp import FpEngine, P, sqrt_p

HERE = os.path.dirname(os.path.abspath(__file__))
NV = 38748


class Fast:
    def __init__(self):
        self.M = Model(); self.E = FpEngine(self.M)
        self.wit = load_assign(os.path.join(HERE, '..', 'best',
                                            'new_instance_partial_39026.json'))
        self.witp = [w % P for w in self.wit]
        self.comp = self.E.comp
        self.v2a = self.E.var2atoms

    def run(self, boolpolicy, forced=None, preassign=None):
        M = self.M; comp = self.comp; v2a = self.v2a
        val = [None] * NV
        reason = [None] * NV
        forced = dict(forced or {})
        for k, x in (preassign or {}).items():
            val[k] = x % P
            reason[k] = 'pre'
        Q = collections.deque(range(M.na))
        inq = bytearray(b'\x01') * M.na
        conflicts = []
        dead = bytearray(M.na)
        decisions = []

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
            val[u] = x; reason[u] = r
            for b in v2a[u]:
                if not inq[b]:
                    inq[b] = 1; Q.append(b)

        while True:
            branch = []
            while Q:
                a = Q.popleft(); inq[a] = 0
                red = reduce_atom(a)
                vs = set()
                for k in red:
                    vs |= set(k)
                if not vs:
                    if red and not dead[a]:
                        dead[a] = 1; conflicts.append(a)
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
                        assign(u, (-c1) * pow(2 * c2 % P, -1, P) % P, a)
                        continue
                    inv = pow(2 * c2 % P, -1, P)
                    roots = sorted({(-c1 + s) * inv % P, (-c1 - s) * inv % P})
                    if len(roots) == 1:
                        assign(u, roots[0], a)
                    else:
                        branch.append((u, a, roots))
            pend = [(u, a, r) for u, a, r in branch if val[u] is None]
            if not pend:
                break
            for u, a, roots in pend:
                if val[u] is not None:
                    continue
                x = boolpolicy(u, roots)
                assign(u, x, 'dec'); decisions.append(u)
        self.reason = reason
        return val, conflicts, decisions

    def eqfail(self, val):
        filled = [0 if x is None else x for x in val]
        E = self.E; M = self.M
        av = [E.eval_atom(a, filled) for a in range(M.na)]
        bad = []
        for e, ts in enumerate(M.eq_terms):
            s = 0
            for c, a in ts:
                s += c * av[a]
            if s % P:
                bad.append(e)
        return bad


if __name__ == '__main__':
    F = Fast()
    pol = sys.argv[1] if len(sys.argv) > 1 else 'wit'
    if pol == 'wit':
        f = lambda u, r: (F.witp[u] if F.witp[u] in r else r[0])
    elif pol == 'zero':
        f = lambda u, r: (0 if 0 in r else r[0])
    else:
        rng = random.Random(int(pol))
        f = lambda u, r: rng.choice(r)
    t = time.time()
    val, conf, dec = F.run(f)
    known = sum(1 for x in val if x is not None)
    print(f"policy={pol} known={known}/{NV} decisions={len(dec)} "
          f"atom-conflicts={len(conf)} t={time.time()-t:.1f}s")
    bad = F.eqfail(val)
    print(f"equations failing mod p (unknowns:=0): {len(bad)}")
    for a in conf[:40]:
        print("   a%d %s" % (a, F.M.src[a][:120]))
    pickle.dump({'val': val, 'conf': conf, 'dec': dec},
                open(os.path.join(HERE, f'bool_{pol}.pkl'), 'wb'))
