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
         chunk=16, dadda_height=2):
    s = p.bit_length()
    Q = MMQB(chunk=chunk, mode=mode, dadda_height=dadda_height)
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


def L1(p, **kw):
    """exhaustive over inputs; returns (#checked, #ok)."""
    square = kw.get('square', False)
    Q, A, B, C, _ = make(p, **kw)
    assert Q.and_weight_ok(), "AND weight does not dominate its local load"
    s = p.bit_length()
    n_ok = n_bad = 0
    pairs = [(a, a) for a in range(p)] if square else \
            [(a, b) for a in range(p) for b in range(p)]
    for a, b in pairs:
        right = a * b % p
        for c in range(1 << s):
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


def L0(p, **kw):
    """total brute force over every variable.  Only for p with <= ~22 vars."""
    square = kw.get('square', False)
    Q, A, B, C, _ = make(p, **kw)
    n = Q.n
    assert n <= 24, f"too many variables for total brute force: {n}"
    s = p.bit_length()
    zeros = set()
    for bitsx in range(1 << n):
        x = [(bitsx >> i) & 1 for i in range(n)]
        if Q.energy(x) == 0:
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


if __name__ == '__main__':
    print("=" * 78)
    print("L0  TOTAL BRUTE FORCE over every variable (ground-state set == truth)")
    print("=" * 78)
    for p in (3, 5, 7):
        for red in ('quotient', 'naf', 'fold'):
            for mode in ('binary', 'wallace', 'dadda', 'unary'):
                try:
                    Q, zeros, truth = L0(p, red=red, mode=mode, leaf=8)
                except AssertionError as ex:
                    print(f"  p={p:2d} {red:9s} {mode:8s}  skipped ({ex})")
                    continue
                ok = (zeros == truth)
                print(f"  p={p:2d} {red:9s} {mode:8s}  vars={Q.n:3d}  "
                      f"|E=0|={len(zeros):3d}  |truth|={len(truth):3d}  "
                      f"{'FAITHFUL' if ok else 'MISMATCH ' + str(zeros ^ truth)}")
    print()
    print("=" * 78)
    print("L1  EXHAUSTIVE OVER INPUTS -- every (a,b,c) on small pseudo-Mersenne p")
    print("=" * 78)
    for p in (13, 29, 61, 127, 251):
        for mult in ('schoolbook', 'karatsuba', 'toom3'):
            for red in ('quotient', 'naf', 'fold'):
                for mode in ('binary', 'wallace'):
                    for sq in (False, True):
                        Q, ok, bad = L1(p, mult=mult, leaf=3, red=red,
                                        mode=mode, square=sq)
                        tag = f"p={p:3d} {mult:10s} {red:9s} {mode:8s} sq={sq:d}"
                        v, cl, jr = Q.triple()
                        print(f"  {tag}  vars={v:6d} K={cl:3d} |J|=2^{jr:<3d} "
                              f"checked={ok + bad:6d} bad={bad}"
                              f"  {'OK' if bad == 0 else '*** FAIL ***'}")
