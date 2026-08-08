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

  3. ROOF DUALITY (QPBO).  Strong persistency labels a variable only if that
     value occurs in a global minimum, so the set roof duality can fix is a
     SUBSET of the ceiling in (2).  The ceiling is therefore the decisive
     measurement and is what is reported; a max-flow QPBO was written and
     discarded after it labelled variables whose labelling had energy 3 (i.e.
     the reachability rule it used was not the sound one), and a wrong
     implementation is worth less than the exact bound it would be bounded by.
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
                lo1 = lo - min(0, c) + c                   # interval with x_v = 1
                hi1 = hi - max(0, c) + c
                lo0 = lo - min(0, c)                       # interval with x_v = 0
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
    print("(d) WHAT the constant variables actually are")
    for p in (29, 61, 127):
        Q, const = ceiling(p, mult='schoolbook', leaf=8, red='naf', mode='wallace')
        kinds = defaultdict(int)
        for v in const:
            kinds[Q.kind[v]] += 1
        print(f"    p={p:4d}: {len(const)} of {Q.n} -- {dict(kinds)}")
    print("""
    They are range slack, nothing else: the top bit of a quotient word or of a
    carry whose declared interval is one wider than its reachable one.  There is
    no propagation into the multiplication itself, and there cannot be: with A
    and B free, every partial product a_i b_j takes both values, so no AND
    ancilla, no carry and no result bit is constant over the ground states.
    Roof duality is bounded by this set, so it cannot do better.

    The factoring literature's table rules (Jiang/Dattani et al.) do fire, but
    only because there the product is a FIXED number and both factors are
    unknown -- the columns are pinned from the top.  Row (c) is that situation
    and it fixes 266 of 200,699 variables: 0.13%.  In the ladder the product of
    each multiplication is another unknown word, so even that does not apply,
    except at the two final comparisons against T.""")


