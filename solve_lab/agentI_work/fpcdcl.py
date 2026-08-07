#!/usr/bin/env python3
"""Mod-p propagation with an implication graph, so conflicts can be traced to
the decisions that caused them (no-good extraction)."""
import pickle, os, collections, sys, random, time
from model import Model, load_assign
from fp import FpEngine, P

HERE = os.path.dirname(os.path.abspath(__file__))
NV = 38748


class Tracer:
    def __init__(self, M, E):
        self.M = M; self.E = E
        self.reason = [None] * NV     # atom index or ('dec',)
        self.level = [-1] * NV

    def analyze(self, atom, val):
        """Return the set of decision variables in the cone of `atom`."""
        seen = set()
        decs = set()
        stack = list(self.E.avarlist[atom])
        while stack:
            v = stack.pop()
            if v in seen or val[v] is None:
                continue
            seen.add(v)
            r = self.reason[v]
            if r is None or r == 'dec':
                decs.add(v)
            else:
                stack.extend(self.E.avarlist[r])
        return decs, seen


def solve(policy='wit', seed=1, max_dec=200000, verbose=True):
    M = Model(); E = FpEngine(M); T = Tracer(M, E)
    wit = load_assign(os.path.join(HERE, '..', 'best',
                                   'new_instance_partial_39026.json'))
    rng = random.Random(seed)
    val = [None] * NV
    unk = [len(E.avarlist[a]) for a in range(M.na)]
    Q = collections.deque(a for a in range(M.na) if unk[a] <= 1)
    decisions = []

    def assign(u, x, reason):
        val[u] = x
        T.reason[u] = reason
        for b in E.var2atoms[u]:
            unk[b] -= 1
            if unk[b] <= 1:
                Q.append(b)

    def propagate():
        branch = []
        while Q:
            a = Q.popleft()
            miss = [x for x in E.avarlist[a] if val[x] is None]
            if len(miss) > 1:
                continue
            if not miss:
                if E.eval_atom(a, val) != 0:
                    return ('conflict', a), branch
                continue
            u = miss[0]
            c2, c1, c0 = E.reduce(a, val, u)
            if c2 == 0:
                if c1 == 0:
                    if c0 != 0:
                        return ('conflict', a), branch
                    continue
                assign(u, (-c0) * pow(c1, -1, P) % P, a)
            else:
                from fp import sqrt_p
                d = (c1 * c1 - 4 * c2 * c0) % P
                s = sqrt_p(d)
                if s is None:
                    return ('conflict', a), branch
                inv = pow(2 * c2 % P, -1, P)
                roots = sorted({(-c1 + s) * inv % P, (-c1 - s) * inv % P})
                if len(roots) == 1:
                    assign(u, roots[0], a)
                else:
                    branch.append((u, a, roots))
        return None, branch

    t0 = time.time()
    while True:
        res, branch = propagate()
        if res:
            _, a = res
            decs, cone = T.analyze(a, val)
            print(f"CONFLICT at a{a}: {M.src[a][:130]}")
            print(f"  cone size {len(cone)}, decisions involved: {len(decs)}")
            print(f"  decision vars: {sorted(decs)[:40]}")
            return val, ('conflict', a, decs, cone), M, E, T
        pend = [(u, aa, r) for u, aa, r in branch if val[u] is None]
        if pend:
            for u, aa, roots in pend:
                if val[u] is not None:
                    continue
                if policy == 'wit':
                    w = wit[u] % P; x = w if w in roots else roots[0]
                elif policy == 'zero':
                    x = 0 if 0 in roots else roots[0]
                else:
                    x = rng.choice(roots)
                assign(u, x, 'dec'); decisions.append(u)
            continue
        unknown = [v for v in range(NV) if val[v] is None]
        if not unknown:
            break
        deg = collections.Counter()
        for a in range(M.na):
            if unk[a] == 2:
                for x in E.avarlist[a]:
                    if val[x] is None:
                        deg[x] += 1
        u = deg.most_common(1)[0][0] if deg else unknown[0]
        if policy == 'wit':
            x = wit[u] % P
        elif policy == 'zero':
            x = 0
        else:
            x = rng.randrange(P)
        assign(u, x, 'dec'); decisions.append(u)
    print(f"COMPLETE mod-p assignment, decisions={len(decisions)}, t={time.time()-t0:.0f}s")
    return val, ('ok', decisions), M, E, T


if __name__ == '__main__':
    pol = sys.argv[1] if len(sys.argv) > 1 else 'wit'
    val, info, M, E, T = solve(pol)
    pickle.dump({'val': val, 'info': info[0]},
                open(os.path.join(HERE, f'fpcdcl_{pol}.pkl'), 'wb'))
