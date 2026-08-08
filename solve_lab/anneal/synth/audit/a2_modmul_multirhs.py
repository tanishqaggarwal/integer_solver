#!/usr/bin/env python3
"""ATTACK 1b/2: faithfulness of build_modmul in the EXACT forms the ladder uses.

verify.py only tests  A*B == C (mod p)  with a single word C.  ladder2.py uses:
    build_modmul(Q,p,lam,d,   [(e,1)])                # 1-word RHS, list form
    build_modmul(Q,p,lam,lam, [(x3,1),(X1,1),(x2,1)])  # SQUARING, 3-word RHS
    build_modmul(Q,p,lam,mm,  [(y3,1),(Y1,1)])         # 2-word RHS
None of the multi-word-RHS / squaring reduction paths are covered by verify.py.
We test them here two ways:
  (L1) exhaustive over all free inputs, canonical witness must have E=0 iff the
       intended congruence holds (and never E=0 when it does not).
  (L0X) full zero-energy ground-state enumeration (validated enumerator) and
       compare its projection to the truth set -- catches spurious completions.
"""
import os, sys, itertools, random
HERE = os.path.dirname(os.path.abspath(__file__))
SQ = os.path.join(HERE, '..', '..', 'squeeze')
sys.path.insert(0, os.path.abspath(SQ))
from mmqb import MMQB
from mm import build_modmul
import verify


def build(p, nrhs, const=0, square=False, mult='schoolbook', red='naf', mode='wallace', leaf=3):
    s = p.bit_length()
    Q = MMQB(chunk=16, mode=mode)
    A = Q.mkword('A', s, lambda wv: wv['_a'])
    B = A if square else Q.mkword('B', s, lambda wv: wv['_b'])
    rhs_words = []
    for i in range(nrhs):
        rhs_words.append(Q.mkword(f'R{i}', s, (lambda wv, i=i: wv[f'_r{i}'])))
    C = [(w, 1) for w in rhs_words]
    build_modmul(Q, p, A, B, C, mult=mult, leaf=leaf, red=red, const=const)
    Q.finalize()
    return Q, A, B, rhs_words


def replay(Q, A, B, rhs_words, a, b, rs, square):
    inp = {}
    for t, v in enumerate(A.bits): inp[v] = (a >> t) & 1
    if not square:
        for t, v in enumerate(B.bits): inp[v] = (b >> t) & 1
    for i, w in enumerate(rhs_words):
        for t, v in enumerate(w.bits): inp[v] = (rs[i] >> t) & 1
    wv0 = {'_a': a, '_b': b}
    for i in range(len(rhs_words)): wv0[f'_r{i}'] = rs[i]
    try:
        x, _ = Q.witness(inp, wv0)
    except (AssertionError, ValueError, ZeroDivisionError):
        return None
    return Q.energy(x)


def L1_multi(p, nrhs, const=0, square=False, sample=None, **kw):
    Q, A, B, rhs_words = build(p, nrhs, const=const, square=square, **kw)
    assert Q.and_weight_ok(), "AND weight fails"
    s = p.bit_length()
    rnd = random.Random(0)
    bad = 0; checked = 0
    apairs = [(a, a) for a in range(p)] if square else [(a, b) for a in range(p) for b in range(p)]
    for a, b in apairs:
        # sample rhs word tuples
        if sample is None:
            rhs_iter = itertools.product(range(1 << s), repeat=nrhs)
        else:
            rhs_iter = [tuple(rnd.randrange(1 << s) for _ in range(nrhs)) for _ in range(sample)]
        for rs in rhs_iter:
            checked += 1
            lhs = a * b
            rhs = sum(rs) + const
            intended = (lhs - rhs) % p == 0
            e = replay(Q, A, B, rhs_words, a, b, rs, square)
            got0 = (e == 0)
            if got0 != intended:
                bad += 1
                if bad <= 5:
                    print(f"    BAD a={a} b={b} rs={rs} const={const} intended={intended} E={e}")
    return Q, checked, bad


def L0X_multi(p, nrhs, const=0, square=False, **kw):
    """full ground-state enumeration; compare projection to truth."""
    Q, A, B, rhs_words = build(p, nrhs, const=const, square=square, **kw)
    s = p.bit_length()
    zeros = set()
    for x in verify.zero_states(Q):
        assert Q.energy(x) == 0
        a = sum(x[v] << t for t, v in enumerate(A.bits))
        b = sum(x[v] << t for t, v in enumerate(B.bits))
        rs = tuple(sum(x[v] << t for t, v in enumerate(w.bits)) for w in rhs_words)
        zeros.add((a, b, rs))
    truth = set()
    rng = range(1 << s)
    for a in rng:
        for b in ([a] if square else rng):
            for rs in itertools.product(rng, repeat=nrhs):
                if (a * b - sum(rs) - const) % p == 0:
                    truth.add((a, b, rs))
    return Q, zeros, truth


if __name__ == '__main__':
    sys.setrecursionlimit(1000000)
    print("=== L1 exhaustive over inputs (multi-word RHS + squaring + const) ===")
    fails = 0
    # small p full exhaustive
    for p in (3, 5, 7):
        for (nrhs, square, const) in [(1,False,0),(2,False,0),(3,True,0),(2,True,0),(1,False,1),(2,False,3)]:
            for red in ('naf','quotient','fold'):
                if red == 'fold' and (nrhs != 1 or const != 0):
                    continue  # fold only supports single-word RHS, const=0 (asserted in mm.py)
                for mode in ('wallace','binary'):
                    sample = None if (nrhs<=2 and p<=5) else 40
                    Q, ch, bad = L1_multi(p, nrhs, const=const, square=square,
                                          sample=sample, red=red, mode=mode)
                    fails += bad
                    print(f"  p={p} nrhs={nrhs} sq={int(square)} const={const} {red:8s} {mode:7s} "
                          f"checked={ch:6d} bad={bad} {'OK' if bad==0 else '*** FAIL ***'}")
    print(f"\n  L1 total bad = {fails}")

    print("\n=== L0X full ground-state enumeration (catches spurious completions) ===")
    fails2 = 0
    for p in (3, 5):
        for (nrhs, square, const) in [(1,False,0),(2,False,0),(3,True,0),(1,False,1)]:
            for red in ('naf','fold'):
                if red == 'fold' and (nrhs != 1 or const != 0):
                    continue
                for mode in ('wallace','binary'):
                    Q, zeros, truth = L0X_multi(p, nrhs, const=const, square=square, red=red, mode=mode)
                    ok = zeros == truth
                    fails2 += 0 if ok else 1
                    extra = sorted(zeros - truth)[:3]
                    print(f"  p={p} nrhs={nrhs} sq={int(square)} const={const} {red:8s} {mode:7s} "
                          f"vars={Q.n:4d} |E=0|={len(zeros):4d} |truth|={len(truth):4d} "
                          f"{'FAITHFUL' if ok else '*** SPURIOUS: '+str(extra)+' ***'}")
    print(f"\n  L0X total mismatches = {fails2}")
    print(f"\nGRAND TOTAL FAILURES: {fails + fails2}")
