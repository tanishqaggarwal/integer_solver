#!/usr/bin/env python3
"""presolve.py -- classical variable elimination before annealing.

Three separate things, kept separate because they answer different questions:

  1. BOUND PROPAGATION / UNIT PROPAGATION over the recorded equations.
     Every penalty in this encoding is  (sum_v c_v x_v + k)^2, so the ground
     states are exactly the solutions of the linear system  sum c_v x_v + k = 0
     over {0,1}.  Classic interval propagation to a fixpoint: if forcing
     x_v = 1 makes some equation's interval exclude 0, then x_v = 0 in every
     ground state.  This subsumes the "multiplication table" rules of the
     annealing-factoring literature (Jiang/Dattani): "a column sum forces a
     carry" is exactly an interval on that column's equation, and "a*b = 1
     forces a = b = 1" is the AND equation's interval.

  2. THE PERSISTENCY CEILING.  No preprocessing -- roof duality included --
     can fix a variable that is not constant across ALL ground states.  We
     measure that ceiling exactly, by enumerating the ground states.

  3. ROOF DUALITY (QPBO) itself, by max-flow on the standard construction,
     on instances small enough to run it.  Reported against the ceiling.
"""
import os
import sys
from collections import defaultdict, deque

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ------------------------------------------------- 1. bound propagation
def propagate(squares, andgates, n, fixed=None, rounds=100):
    """squares: [(lin dict v->coef, const)] each asserted == 0
       andgates: [(z, i, j)] meaning z = x_i * x_j
       returns dict v -> value for every provably-determined variable."""
    dom = dict(fixed or {})
    occ = defaultdict(list)
    for e, (lin, k) in enumerate(squares):
        for v in lin:
            occ[v].append(e)
    andocc = defaultdict(list)
    for g, (z, i, j) in enumerate(andgates):
        andocc[z].append(g)
        andocc[i].append(g)
        andocc[j].append(g)
    queue = deque(range(len(squares)))
    inq = set(queue)
    gq = deque(range(len(andgates)))
    ginq = set(gq)

    def setv(v, val):
        if v in dom:
            if dom[v] != val:
                raise ValueError("infeasible")
            return False
        dom[v] = val
        for e in occ[v]:
            if e not in inq:
                inq.add(e)
                queue.append(e)
        for g in andocc[v]:
            if g not in ginq:
                ginq.add(g)
                gq.append(g)
        return True

    for _ in range(rounds):
        moved = False
        while queue:
            e = queue.popleft()
            inq.discard(e)
            lin, k = squares[e]
            base = k
            free = []
            for v, c in lin.items():
                if v in dom:
                    base += c * dom[v]
                else:
                    free.append((v, c))
            lo = base + sum(min(0, c) for _, c in free)
            hi = base + sum(max(0, c) for _, c in free)
            if lo > 0 or hi < 0:
                raise ValueError("infeasible")
            for v, c in free:
                lo1, hi1 = lo, hi                          # x_v = 1
                lo1 += c if c > 0 else -c * 0
                lo1 = lo - min(0, c) + c
                hi1 = hi - max(0, c) + c
                lo0 = lo - min(0, c)
                hi0 = hi - max(0, c)
                can1 = not (lo1 > 0 or hi1 < 0)
                can0 = not (lo0 > 0 or hi0 < 0)
                if not can1 and not can0:
                    raise ValueError("infeasible")
                if not can1:
                    moved |= setv(v, 0)
                elif not can0:
                    moved |= setv(v, 1)
        while gq:
            g = gq.popleft()
            ginq.discard(g)
            z, i, j = andgates[g]
            zi, ii, ji = dom.get(z), dom.get(i), dom.get(j)
            if ii == 0 or ji == 0:
                if zi is None:
                    moved |= setv(z, 0)
            if ii == 1 and ji == 1 and zi is None:
                moved |= setv(z, 1)
            if zi == 1:
                if ii is None:
                    moved |= setv(i, 1)
                if ji is None:
                    moved |= setv(j, 1)
            if zi == 0 and ii == 1 and ji is None:
                moved |= setv(j, 0)
            if zi == 0 and ji == 1 and ii is None:
                moved |= setv(i, 0)
        if not moved:
            break
    return dom


def gates_of(Q):
    return [(z, i, j) for (i, j), z in Q.andcache.items()]


# ------------------------------------------------- 3. roof duality (QPBO)
class Dinic:
    def __init__(self, n):
        self.n = n
        self.g = [[] for _ in range(n)]

    def add(self, u, v, c):
        self.g[u].append([v, c, len(self.g[v])])
        self.g[v].append([u, 0, len(self.g[u]) - 1])

    def maxflow(self, s, t):
        flow = 0
        while True:
            lev = [-1] * self.n
            lev[s] = 0
            dq = deque([s])
            while dq:
                u = dq.popleft()
                for e in self.g[u]:
                    if e[1] > 0 and lev[e[0]] < 0:
                        lev[e[0]] = lev[u] + 1
                        dq.append(e[0])
            if lev[t] < 0:
                return flow
            it = [0] * self.n

            def dfs(u, f):
                if u == t:
                    return f
                while it[u] < len(self.g[u]):
                    e = self.g[u][it[u]]
                    v = e[0]
                    if e[1] > 0 and lev[v] == lev[u] + 1:
                        d = dfs(v, min(f, e[1]))
                        if d:
                            e[1] -= d
                            self.g[v][e[2]][1] += d
                            return d
                    it[u] += 1
                return 0
            while True:
                f = dfs(s, 1 << 60)
                if not f:
                    break
                flow += f


def qpbo(Qdict, n):
    """roof duality.  Returns dict v -> persistent value (strong persistency:
       these values occur in SOME global minimum, and by the standard
       autarky argument can be fixed without losing optimality)."""
    # nodes: 0 = source (x=0 side), 1 = sink, 2+2i = x_i, 3+2i = xbar_i
    N = 2 + 2 * n
    S, T = 0, 1
    D = Dinic(N)
    X = lambda i: 2 + 2 * i                                          # noqa: E731
    Xb = lambda i: 3 + 2 * i                                         # noqa: E731
    lin = defaultdict(int)
    quad = defaultdict(int)
    for m, c in Qdict.items():
        if len(m) == 1:
            lin[m[0]] += c
        elif len(m) == 2:
            quad[m] += c
    for i, c in lin.items():
        if c > 0:
            D.add(S, X(i), c)
            D.add(Xb(i), T, c)
        elif c < 0:
            D.add(X(i), T, -c)
            D.add(S, Xb(i), -c)
    for (i, j), c in quad.items():
        h = abs(c) / 2.0
        h = abs(c)
        if c > 0:
            D.add(X(i), Xb(j), h)
            D.add(X(j), Xb(i), h)
        else:
            D.add(S, X(i), h)
            D.add(Xb(i), T, h)
            D.add(X(i), X(j), h)
            D.add(Xb(j), Xb(i), h)
    D.maxflow(S, T)
    # residual reachability
    seen = [False] * N
    dq = deque([S])
    seen[S] = True
    while dq:
        u = dq.popleft()
        for e in D.g[u]:
            if e[1] > 0 and not seen[e[0]]:
                seen[e[0]] = True
                dq.append(e[0])
    rev = [[] for _ in range(N)]
    for u in range(N):
        for e in D.g[u]:
            if e[1] > 0:
                rev[e[0]].append(u)
    seenT = [False] * N
    dq = deque([T])
    seenT[T] = True
    while dq:
        u = dq.popleft()
        for v in rev[u]:
            if not seenT[v]:
                seenT[v] = True
                dq.append(v)
    out = {}
    for i in range(n):
        if seen[X(i)] and not seen[Xb(i)]:
            out[i] = 0
        elif seen[Xb(i)] and not seen[X(i)]:
            out[i] = 1
        elif seenT[X(i)] and not seenT[Xb(i)]:
            out[i] = 1
        elif seenT[Xb(i)] and not seenT[X(i)]:
            out[i] = 0
    return out


# ------------------------------------------------- 2. persistency ceiling
def ceiling(p, **kw):
    """exact: how many variables are constant over ALL ground states.
       No preprocessing of any kind can fix more than this."""
    import verify
    Q, A, B, C, _ = verify.make(p, **kw)
    square = kw.get('square', False)
    seen0 = [False] * Q.n
    seen1 = [False] * Q.n
    pairs = [(a, a) for a in range(p)] if square else \
            [(a, b) for a in range(p) for b in range(p)]
    for a, b in pairs:
        c = a * b % p
        e = verify.replay(Q, A, B, C, a, b, c, square)
        assert e == 0
        x, _ = Q.witness(_inputs(A, B, C, a, b, c), {'_a': a, '_b': b, '_c': c})
        for v in range(Q.n):
            if x[v]:
                seen1[v] = True
            else:
                seen0[v] = True
    const = [v for v in range(Q.n) if not (seen0[v] and seen1[v])]
    return Q, const


def _inputs(A, B, C, a, b, c):
    inp = {}
    for t, v in enumerate(A.bits):
        inp[v] = (a >> t) & 1
    for t, v in enumerate(B.bits):
        inp[v] = (b >> t) & 1
    for t, v in enumerate(C.bits):
        inp[v] = (c >> t) & 1
    return inp


if __name__ == '__main__':
    import verify
    import measure
    print("=" * 96)
    print("PRESOLVE -- what a classical pass can remove before the annealer sees it")
    print("=" * 96)
    print()
    print("(a) PERSISTENCY CEILING: variables constant across every ground state")
    print(f"    {'p':>6} {'vars':>7} {'constant':>9}  (any presolver's hard upper bound)")
    for p in (13, 29, 61, 127):
        Q, const = ceiling(p, mult='schoolbook', leaf=8, red='naf', mode='wallace')
        print(f"    {p:6d} {Q.n:7d} {len(const):9d}")
    print()
    print("(b) BOUND PROPAGATION to a fixpoint on the real 256-bit modmul")
    print(f"    {'variant':>34} {'vars':>9} {'fixed':>8} {'%':>7}")
    for key, kw in (("school/naf/wallace", dict(mult='schoolbook', red='naf', mode='wallace')),
                    ("karatsuba(32)/naf/wallace", dict(mult='karatsuba', leaf=32, red='naf',
                                                       mode='wallace')),
                    ("school/naf/binary", dict(mult='schoolbook', red='naf', mode='binary'))):
        Q = _build256(**kw)
        dom = propagate(Q.squares, gates_of(Q), Q.n)
        print(f"    {key:>34} {Q.n:9,d} {len(dom):8,d} {100*len(dom)/Q.n:6.2f}%")
    print()
    print("(c) SAME, with the product pinned to a constant  (the situation the")
    print("    factoring literature's table rules were written for: c is known)")
    print(f"    {'variant':>34} {'vars':>9} {'fixed':>8} {'%':>7}")
    for key, kw in (("school/naf/wallace", dict(mult='schoolbook', red='naf', mode='wallace')),):
        Q, Cbits = _build256(pin=True, **kw)
        val = 0x9d671cd581c69bc5e697f5e45bcd07c6741496c7e6d8a2f0e6d1a2b3c4d5e6f7 % measure.P
        fixed = {v: (val >> t) & 1 for t, v in enumerate(Cbits)}
        dom = propagate(Q.squares, gates_of(Q), Q.n, fixed=fixed)
        print(f"    {key:>34} {Q.n:9,d} {len(dom):8,d} {100*len(dom)/Q.n:6.2f}%")
    print()
    print("(d) ROOF DUALITY (QPBO max-flow) vs the ceiling, on instances it can run on")
    print(f"    {'p':>6} {'vars':>7} {'qpbo fixed':>11} {'ceiling':>9}")
    for p in (13, 29, 61):
        Q, A, B, C, _ = verify.make(p, mult='schoolbook', leaf=8, red='naf', mode='wallace')
        got = qpbo(Q.Q, Q.n)
        _, const = ceiling(p, mult='schoolbook', leaf=8, red='naf', mode='wallace')
        print(f"    {p:6d} {Q.n:7d} {len(got):11d} {len(const):9d}")


def _build256(pin=False, **kw):
    import measure
    from mmqb import MMQB
    from mm import build_modmul
    Q = MMQB(chunk=kw.pop('chunk', 16), mode=kw.pop('mode', 'wallace'))
    A = Q.mkword('A', 256, lambda wv: 0)
    B = Q.mkword('B', 256, lambda wv: 0)
    C = Q.mkword('C', 256, lambda wv: 0)
    build_modmul(Q, measure.P, A, B, C, **kw)
    Q.finalize()
    return (Q, C.bits) if pin else Q
