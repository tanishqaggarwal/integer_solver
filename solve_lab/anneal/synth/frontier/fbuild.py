#!/usr/bin/env python3
"""fbuild.py -- builders for the frontier study.  Imports the existing model /
squeeze encoders unmodified and exposes every modmul instance as a numpy Ising
(via model.qubo_to_ising) plus the metadata a solver harness needs.

Two encoders for ONE modular multiply  a*b == c (mod p):
  * baseline  : model.build_modmul  (ladder.py / qubo.py: wallace|binary, W_and)
  * squeeze   : squeeze/mm.build_modmul on an MMQB (karatsuba+NAF, dadda, ...)

Plus a one-operand-known variant: clamp the bits of `a` to the planted value.
"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'synth', 'solver'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'squeeze'))
import random
import numpy as np
import model as M


# ----------------------------------------------------------- prime picker
def pick_prime(s, seed=3):
    """same recipe model.build_modmul uses, exposed so both encoders share p,a,b,c."""
    rnd = random.Random(seed)
    p = (1 << (s - 1)) + 2 * rnd.randrange(1 << (s - 3)) + 1
    while True:
        if all(p % q for q in range(3, int(p ** .5) + 1, 2)):
            break
        p += 2
    a, b = rnd.randrange(2, p), rnd.randrange(2, p)
    cc = a * b % p
    return p, a, b, cc


# ----------------------------------------------------------- baseline encoder
def baseline_modmul(s, mode='wallace', W_and=None, seed=3):
    mm = M.build_modmul(s, mode=mode, seed=seed, W_and=W_and)
    mm['encoder'] = f"baseline/{mode}" + (f"/Wand{W_and}" if W_and else "")
    return mm


# ----------------------------------------------------------- squeeze encoder
def squeeze_modmul(s, mult='karatsuba', red='naf', leaf=8, mode='wallace', seed=3):
    """Build a*b==c (mod p) with the squeeze encoder (squeeze/mm.py) and return a
    dict shaped like model.build_modmul's output (Q, ising, A, B, p, a, b, c, ...)."""
    from mmqb import MMQB
    import mm as SQ
    p, a, b, cc = pick_prime(s, seed=seed)
    Q = MMQB(mode=mode)
    A = Q.mkword("a", s, lambda wv, a=a: a)
    Bw = Q.mkword("b", s, lambda wv, b=b: b)
    C = Q.mkword("c", s, lambda wv, cc=cc: cc)
    SQ.build_modmul(Q, p, A, Bw, C, mult=mult, leaf=leaf, red=red, tag="mm")
    Q.finalize()
    x, _ = Q.witness({}, {})
    e0 = Q.energy(x)
    assert e0 == 0, f"squeeze planted state energy {e0} != 0 (s={s})"
    ising = M.qubo_to_ising(Q.Q, Q.n)
    sstar = 2 * np.array(x, dtype=np.float64) - 1
    assert abs(ising.energy(sstar)) < 1e-6
    return dict(Q=Q, ising=ising, A=A.bits, B=Bw.bits, Cbits=C.bits,
                p=p, a=a, b=b, c=cc, xstar=x, s=s,
                encoder=f"squeeze/{mult}/{red}/{mode}")


# ----------------------------------------------------------- clamp one operand
def clamp_operand(mm, which='a', rng=None):
    """Return (x0, clamp) that fixes the bits of operand `which` (a or b) to its
    planted value; every other bit random.  The solver then searches only the
    other operand and the ancillas -- the 'one operand known' sub-instance."""
    Q = mm['Q']
    n = Q.n
    rng = rng or np.random.default_rng(0)
    x0 = list(rng.integers(0, 2, size=n).astype(int))
    bits = mm['A'] if which == 'a' else mm['B']
    val = mm['a'] if which == 'a' else mm['b']
    for t, v in enumerate(bits):
        x0[v] = (val >> t) & 1
    clamp = set(bits)
    return x0, clamp


if __name__ == '__main__':
    # smoke test: both encoders produce a valid E=0 planted ground state
    print("s   baseline_n  squeeze_kara_n  squeeze_dadda_n")
    for s in range(4, 13):
        bm = baseline_modmul(s)
        sk = squeeze_modmul(s, mult='karatsuba', red='naf', leaf=8, mode='wallace')
        sd = squeeze_modmul(s, mult='karatsuba', red='naf', leaf=8, mode='dadda')
        print(f"{s:2d}  {bm['Q'].n:8d}   {sk['Q'].n:10d}   {sd['Q'].n:10d}")
    print("OK: all planted states are E=0")
