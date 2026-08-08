#!/usr/bin/env python3
"""common.py -- shared builders for the convex/submodular study.

Everything here IMPORTS the read-only encoders (squeeze/mm.py, squeeze/mmqb.py,
squeeze/verify.py, ../solver/model.py); nothing here modifies them.

Two modmul flavours are used:

  * `build_mm(s, p=None)` -- squeeze/mm.py schoolbook modmul  A*B == C (mod p)
    with C a FREE word (the operands AND the product are unknown).  This is the
    natural "nothing pinned" instance and the one the persistency ceiling in
    squeeze/presolve.py was measured on.  Works at any size because mm.py takes
    the modulus as an argument (model.build_modmul trial-divides and dies > 40b).

  * verify.make(p, ...) -- the same encoder, used by the exhaustive
    ground-state enumerator verify.zero_states for the soundness gate.

A QUBO here is the finalized dict Q.Q : monomial(tuple) -> int coefficient.
Submodular coupler  <=>  quadratic coefficient <= 0  (ferromagnetic / min-cut
minimisable);  supermodular  <=>  coefficient > 0.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SQUEEZE = os.path.abspath(os.path.join(HERE, '..', '..', 'squeeze'))
ANNEAL = os.path.abspath(os.path.join(HERE, '..', '..'))
SOLVER = os.path.abspath(os.path.join(HERE, '..', 'solver'))
for p in (SQUEEZE, ANNEAL, SOLVER):
    if p not in sys.path:
        sys.path.insert(0, p)


def _is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_of_bits(s, seed=1):
    """a prime with exactly s bits (top bit set)."""
    if s <= 2:
        return {1: 1, 2: 3}.get(s, 3)
    import random
    rnd = random.Random(seed * 1000 + s)
    base = 1 << (s - 1)
    while True:
        cand = base + rnd.randrange(1 << (s - 1)) | 1
        cand |= (1 << (s - 1))
        if cand.bit_length() != s:
            continue
        if _is_prime(cand):
            return cand


# the real curve prime (squeeze/measure.py) for the s=256 headline
def real_p256():
    return 2 ** 256 - 2 ** 32 - 977


def build_mm(s, p=None, mult='schoolbook', leaf=32, red='naf', mode='wallace',
             free_c=True):
    """squeeze/mm.py modmul at operand width s.  Returns dict with Q, A, B, C."""
    from mmqb import MMQB
    import mm
    if p is None:
        p = prime_of_bits(s)
    Q = MMQB(chunk=16, mode=mode)
    A = Q.mkword('A', s, lambda wv: wv.get('_a', 0))
    B = Q.mkword('B', s, lambda wv: wv.get('_b', 0))
    if free_c:
        C = Q.mkword('C', s, lambda wv: wv.get('_c', 0))
        mm.build_modmul(Q, p, A, B, C, mult=mult, leaf=leaf, red=red)
    else:
        C = None
        raise NotImplementedError
    Q.finalize()
    return dict(Q=Q, A=A, B=B, C=C, p=p, s=s)


# ------------------------------------------------------------- classification
def coupler_split(Q):
    """counts of submodular (coef<=0) vs supermodular (coef>0) couplers, plus a
    kind x kind breakdown of the supermodular ones.  Uses Q.Q (finalized)."""
    from collections import Counter
    sub = sup = 0
    sup_kinds = Counter()
    sub_kinds = Counter()
    for m, c in Q.Q.items():
        if len(m) != 2:
            continue
        ka, kb = Q.kind[m[0]], Q.kind[m[1]]
        key = tuple(sorted((ka, kb)))
        if c > 0:
            sup += 1
            sup_kinds[key] += 1
        else:
            sub += 1
            sub_kinds[key] += 1
    return dict(sub=sub, sup=sup, total=sub + sup,
                sup_frac_of_super=sup, sub_kinds=dict(sub_kinds),
                sup_kinds=dict(sup_kinds))
