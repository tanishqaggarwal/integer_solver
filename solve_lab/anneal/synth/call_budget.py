#!/usr/bin/env python3
"""call_budget.py -- what does a budget of 2^32 real annealer calls actually buy?

Two independent facts, both measured elsewhere in this repo, combine into a law:

  (1) interval split: outer runs = 2^(b-mu), each call a mu-bit sub-instance.
  (2) NO GRADIENT (synth/solver): one anneal of a sub-instance finds its answer
      only by chance, so anneals-per-sub-instance ~ 2^mu (worse, empirically).
      => total real-annealer calls to recover a b-bit key ~ 2^b, independent of mu.

  (3) EXECUTABILITY: every comb window does full s-bit FIELD arithmetic (s = curve
      bit-size) no matter how few scalar bits mu it resolves. And field arithmetic
      is itself gradient-free -- a*b=c is a needle just like k*G=T. So a call only
      succeeds if the field is inside the 'arithmetic-annealing frontier' F.

Consequence: a budget of N=2^32 calls covers a b-bit SCALAR search only if b<=32,
AND only if the s=b-bit field arithmetic anneals (s <= F). If F < 32 the 2^32
calls cannot be executed successfully at all -- they fail on the arithmetic.
"""
import math

def total_calls(b, gradient=False, split='interval'):
    """real annealer calls to recover a b-bit key."""
    if gradient:
        # a gradient-ful sub-problem would let one call solve mu bits; then the
        # optimal split is mu as large as fits, and total = 2^(b-mu). The ONLY
        # known gradient-ful reformulation is modular subset-sum from relations,
        # which needs index calculus -- absent for prime-field ECC. See §subset_sum.
        return None
    return 2 ** b            # no gradient: 2^(b-mu) * 2^mu

if __name__ == '__main__':
    print("Real-annealer calls to recover a b-bit key (gradient-free encoding):")
    print(f"  {'b':>4} {'calls=2^b':>10} {'<= 2^32 budget?':>16}")
    for b in (8, 16, 24, 32, 40, 48, 64, 128, 256):
        print(f"  {b:4d} {'2^%d'%b:>10} {'YES' if b<=32 else 'no':>16}")
    print("""
So in the COUNTING sense 2^32 calls reaches a 32-bit key -- but only if the
32-bit field arithmetic anneals.  The binding question is the arithmetic
frontier F: the largest field size whose modular multiply an annealer actually
settles.  If F ~ 8, no 32-bit instance is executable and 2^32 buys nothing above
an 8-bit key (which needs only ~2^8 calls).  Measuring F is the whole game.

The one escape is a GRADIENT-FUL sub-problem, where one call solves many bits and
total drops to ~2^(b-mu). The only known such reformulation is modular
subset-sum from a factor base of relations (energy (sum c_j l_j - target)^2 is
LINEAR in the bits => real gradient). It needs relations with known logs, i.e.
index calculus -- which does not exist for prime-field ECC. That absence is
exactly why prime-field ECDLP is classically hard, and it is what denies the
annealer a gradient here.
""")
