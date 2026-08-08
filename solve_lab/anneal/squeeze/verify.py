#!/usr/bin/env python3
"""verify.py -- faithfulness of every modmul variant, exhaustively.

Three layers, strongest first:

  L0  TOTAL brute force.  Tiny p (2 or 3 bits): enumerate EVERY assignment of
      EVERY variable in the Hamiltonian and check
          { x : E(x) = 0 }  ==  { x : a*b == c (mod p) and all ancillas correct }.
      No structural argument, no replay -- the actual ground-state set.

  L1  Exhaustive over inputs.  For every (a,b,c) in [0,p)x[0,p)x[0,2^s):
        * if c == a*b mod p the replayed witness must have E = 0
        * otherwise no witness exists (replay raises) AND, because every
          penalty is a perfect square of an integer form and the AND weight
          dominates its own local load (checked by and_weight_ok), no
          completion can reach 0 either.  We test that claim directly at L0.

  L2  The whole-instance test the lab already uses: rebuild the windowed comb
      ladder with the new modmul and re-run demo_win's exhaustive scan over
      every candidate scalar on a small curve (see demo_win2.py).
"""
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmqb import MMQB                                            # noqa: E402
from mm import build_modmul                                      # noqa: E402


def make(p, mult='schoolbook', leaf=8, red='naf', mode='wallace', square=False,
         chunk=16, dadda_height=2, naf_merge=True):
    s = p.bit_length()
    Q = MMQB(chunk=chunk, mode=mode, dadda_height=dadda_height,
             naf_merge=naf_merge)
    A = Q.mkword('A', s, lambda wv: wv['_a'])
    B = A if square else Q.mkword('B', s, lambda wv: wv['_b'])
    C = Q.mkword('C', s, lambda wv: wv['_c'])
    base = Q.n
    build_modmul(Q, p, A, B, C, mult=mult, leaf=leaf, red=red)
    Q.finalize()
    return Q, A, B, C, base


def replay(Q, A, B, C, a, b, c, square):
    inp = {}
    for t, v in enumerate(A.bits):
        inp[v] = (a >> t) & 1
    for t, v in enumerate(B.bits):
        inp[v] = (b >> t) & 1
    for t, v in enumerate(C.bits):
        inp[v] = (c >> t) & 1
    wv0 = {'_a': a, '_b': b, '_c': c}
    try:
        x, _ = Q.witness(inp, wv0)
    except AssertionError:
        return None
    return Q.energy(x)


def L1(p, sample=None, **kw):
    """exhaustive over inputs; returns (#checked, #ok).
       sample=None scans every c in [0,2^s); an int scans the correct c plus
       that many random wrong ones (for p where the full scan is too slow)."""
    import random
    square = kw.get('square', False)
    Q, A, B, C, _ = make(p, **kw)
    assert Q.and_weight_ok(), "AND weight does not dominate its local load"
    s = p.bit_length()
    n_ok = n_bad = 0
    rnd = random.Random(0)
    pairs = [(a, a) for a in range(p)] if square else \
            [(a, b) for a in range(p) for b in range(p)]
    for a, b in pairs:
        right = a * b % p
        cs = range(1 << s) if sample is None else \
            [right] + [rnd.randrange(1 << s) for _ in range(sample)]
        for c in cs:
            e = replay(Q, A, B, C, a, b, c, square)
            if c == right:
                if e != 0:
                    n_bad += 1
                else:
                    n_ok += 1
            else:
                # a non-canonical duplicate c = right + p < 2^s is a legal
                # zero-energy state of the base encoding too; count it as such.
                if e == 0 and c != right + p:
                    n_bad += 1
                else:
                    n_ok += 1
    return Q, n_ok, n_bad


def zero_states(Q):
    """Enumerate EVERY assignment with E = 0, exactly, without touching 2^n.

    E is a sum of squares of integer linear forms, of Rosenberg AND penalties,
    and (unary mode only) of thermometer terms  c_{t+1}(1-c_t) -- all
    non-negative.  So  E(x) = 0  iff every one of those is individually zero,
    i.e. x solves the constraint system.  Depth-first search with interval
    propagation enumerates that solution set in time proportional to its size,
    which is what makes an exhaustive statement possible past 22 variables."""
    eqs = Q.squares
    gates = [(z, i, j) for (i, j), z in Q.andcache.items()]
    orders = Q.orders
    n = Q.n
    occ = [[] for _ in range(n)]
    for e, (lin, k) in enumerate(eqs):
        for v in lin:
            occ[v].append(e)
    out = []
    val = [None] * n

    def feasible(elist):
        for e in elist:
            lin, k = eqs[e]
            lo = hi = k
            for v, c in lin.items():
                if val[v] is None:
                    lo += min(0, c)
                    hi += max(0, c)
                else:
                    lo += c * val[v]
                    hi += c * val[v]
            if lo > 0 or hi < 0:
                return False
        return True

    def rec(v):
        if v == n:
            out.append(list(val))
            return
        for b in (0, 1):
            val[v] = b
            ok = feasible(occ[v])
            if ok:
                for (z, i, j) in gates:
                    a, bb, cc = val[i], val[j], val[z]
                    if a is not None and bb is not None and cc is not None and cc != a * bb:
                        ok = False
                        break
                    if cc == 1 and ((a == 0) or (bb == 0)):
                        ok = False
                        break
            if ok:
                for (u, w) in orders:
                    if val[u] == 0 and val[w] == 1:
                        ok = False
                        break
            if ok:
                rec(v + 1)
        val[v] = None
    rec(0)
    return out


def L0X(p, **kw):
    """exhaustive ground-state enumeration via zero_states (no 2^n scan)."""
    square = kw.get('square', False)
    Q, A, B, C, _ = make(p, **kw)
    s = p.bit_length()
    zeros = set()
    for x in zero_states(Q):
        assert Q.energy(x) == 0
        a = sum(x[v] << t for t, v in enumerate(A.bits))
        b = sum(x[v] << t for t, v in enumerate(B.bits))
        c = sum(x[v] << t for t, v in enumerate(C.bits))
        zeros.add((a, b, c))
    truth = set()
    for a in range(1 << s):
        for b in range(1 << s):
            if square and a != b:
                continue
            for c in range(1 << s):
                if (a * b - c) % p == 0:
                    truth.add((a, b, c))
    return Q, zeros, truth


def L0(p, cap=20, **kw):
    """total brute force over every variable.  Only for tiny p."""
    square = kw.get('square', False)
    Q, A, B, C, _ = make(p, **kw)
    n = Q.n
    assert n <= cap, f"{n} variables > cap {cap}"
    s = p.bit_length()
    off = Q.Q.get((), 0)
    lin = [(m[0], c) for m, c in Q.Q.items() if len(m) == 1]
    qua = [(m[0], m[1], c) for m, c in Q.Q.items() if len(m) == 2]
    zeros = set()
    for bitsx in range(1 << n):
        e = off
        for v, c in lin:
            if (bitsx >> v) & 1:
                e += c
        for u, v, c in qua:
            if (bitsx >> u) & 1 and (bitsx >> v) & 1:
                e += c
        if e == 0:
            a = sum(((bitsx >> v) & 1) << t for t, v in enumerate(A.bits))
            b = sum(((bitsx >> v) & 1) << t for t, v in enumerate(B.bits))
            c = sum(((bitsx >> v) & 1) << t for t, v in enumerate(C.bits))
            zeros.add((a, b, c))
    truth = set()
    for a in range(1 << s):
        for b in range(1 << s):
            if square and a != b:
                continue
            for c in range(1 << s):
                if (a * b - c) % p == 0:
                    truth.add((a, b, c))
    return Q, zeros, truth


if __name__ == '__main__':
    import sys as _s
    _s.setrecursionlimit(100000)
    fails = 0
    print("=" * 92)
    print("L0  EXHAUSTIVE GROUND STATES -- every zero-energy assignment enumerated,")
    print("    compared against the truth set { (a,b,c) : a*b == c (mod p) }.")
    print("=" * 92)
    for p in (3, 5, 7):
        for mult in ('schoolbook', 'karatsuba', 'toom3'):
            for red in ('quotient', 'naf', 'fold'):
                for mode in ('binary', 'wallace', 'dadda', 'unary'):
                    Q, zeros, truth = L0X(p, mult=mult, leaf=3, red=red, mode=mode)
                    ok = (zeros == truth)
                    fails += 0 if ok else 1
                    print(f"  p={p:2d} {mult:10s} {red:9s} {mode:8s} vars={Q.n:4d}"
                          f"  |E=0|={len(zeros):4d}  |truth|={len(truth):4d}  "
                          f"{'FAITHFUL' if ok else 'MISMATCH ' + str(sorted(zeros ^ truth))[:90]}")
    print()
    print("  and the same at p=13 (4-bit words), the widest the full enumeration reaches:")
    for mult in ('schoolbook', 'karatsuba', 'toom3'):
        for mode in ('binary', 'wallace'):
            Q, zeros, truth = L0X(13, mult=mult, leaf=3, red='naf', mode=mode)
            ok = (zeros == truth)
            fails += 0 if ok else 1
            print(f"  p=13 {mult:10s} naf       {mode:8s} vars={Q.n:4d}"
                  f"  |E=0|={len(zeros):4d}  |truth|={len(truth):4d}  "
                  f"{'FAITHFUL' if ok else 'MISMATCH'}")
    print()
    print("=" * 92)
    print("L1  EXHAUSTIVE OVER INPUTS -- every (a,b,c) triple, every variant")
    print("=" * 92)
    for p in (13, 29, 61):
        for mult in ('schoolbook', 'karatsuba', 'toom3'):
            for red in ('quotient', 'naf', 'fold'):
                for mode in ('binary', 'wallace', 'dadda', 'unary'):
                    for sq in (False, True):
                        Q, ok, bad = L1(p, mult=mult, leaf=3, red=red,
                                        mode=mode, square=sq)
                        fails += bad
                        v, cl, jr = Q.triple()
                        print(f"  p={p:3d} {mult:10s} {red:9s} {mode:8s} sq={sq:d}"
                              f"  vars={v:6d} K={cl:3d} |J|=2^{jr:<3d} "
                              f"checked={ok + bad:7d} bad={bad}"
                              f"  {'OK' if bad == 0 else '*** FAIL ***'}")
    print()
    print("  larger p (all (a,b); correct c plus 16 random wrong c each):")
    for p in (127, 251, 1021, 8191):
        for mult in ('schoolbook', 'karatsuba', 'toom3', [('toom3', 8), ('karatsuba', 4)]):
            for mode in ('binary', 'wallace'):
                nm = mult if isinstance(mult, str) else 'toom3+kara'
                Q, ok, bad = L1(p, sample=16, mult=mult, leaf=4, red='naf', mode=mode)
                fails += bad
                v, cl, jr = Q.triple()
                print(f"  p={p:5d} {nm:10s} naf       {mode:8s}"
                      f"  vars={v:6d} K={cl:3d} |J|=2^{jr:<3d} "
                      f"checked={ok + bad:7d} bad={bad}"
                      f"  {'OK' if bad == 0 else '*** FAIL ***'}")
    print()
    print(f"TOTAL FAILURES: {fails}")
