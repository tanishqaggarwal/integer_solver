#!/usr/bin/env python3
"""lazy.py -- what deferred (lazy) modular reduction actually costs here.

Two numbers decide it, and both are measured, not modelled:

  (1) how much of a modular multiplication IS the reduction, i.e. the price of
      the explicit quotient word -- that is the most lazy reduction can ever
      save;
  (2) how fast a multiplication grows when its operands are allowed to grow --
      that is what lazy reduction has to pay.

In this circuit every deferred value is immediately fed into another
multiplication (x3 -> d -> lam*d), so widths compose MULTIPLICATIVELY: skipping
j reductions makes the operands ~2^j times wider, not j times wider.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import random                              # noqa: E402
from enc import Ladder2                    # noqa: E402
from marginal import modmul, _prime        # noqa: E402


def unreduced_mul(s, wa, wb, mode='binary', kdepth=0, kmin=8, chunk=16, seed=1):
    """A*B == Wfull  exactly over Z: no modulus, no quotient word."""
    rnd = random.Random(seed)
    p = _prime(s, rnd)
    L = Ladder2(p, chunk=chunk, mode=mode, kdepth=kdepth, kmin=kmin)
    Q = L.qb
    zero = lambda wv: 0
    A = Q.word("A", wa, zero)
    Bw = Q.word("B", wb, zero)
    L._prod_word(A, zero, Bw, zero, "u", kdepth)
    Q.finalize()
    st = Q.stats()
    st['vars'] -= (wa + wb)
    return st


if __name__ == '__main__':
    s = 256
    print("(1) how much of a modular multiplication is the REDUCTION")
    print(f"{'variant':>44} {'vars':>10} {'couplers':>13}")
    for kd in (0, 4):
        r = modmul(s, kdepth=kd)
        u = unreduced_mul(s, s, s, kdepth=kd)
        print(f"{'A*B == W (mod p)      kdepth=%d' % kd:>44} "
              f"{r['vars']:10,d} {r['couplers']:13,d}")
        print(f"{'A*B == Wfull over Z   kdepth=%d' % kd:>44} "
              f"{u['vars']:10,d} {u['couplers']:13,d}")
        print(f"{'--> reduction glue is':>44} "
              f"{r['vars'] - u['vars']:+10,d} {r['couplers'] - u['couplers']:+13,d}"
              f"   ({100.0 * (r['vars'] - u['vars']) / r['vars']:.1f}% of the multiply)")
    print()
    print("(2) what wider operands cost (this is the bill lazy reduction pays)")
    print(f"{'|A| x |B|':>16} {'vars':>10} {'couplers':>13} {'AND':>10} {'vs 256x256':>11}")
    base = None
    for wa in (256, 264, 288, 320, 384, 512):
        st = modmul(s, kdepth=0, wa=wa, wb=wa)
        base = base or st['vars']
        print(f"{'%dx%d' % (wa, wa):>16} {st['vars']:10,d} {st['couplers']:13,d} "
              f"{st['and_vars']:10,d} {st['vars'] / base:10.2f}x", flush=True)
    print()
    print("(3) the realistic one-step-deferred shape: x3 left at 512 bits, so the")
    print("    next window's d is 513 bits and lam*d is 256 x 513")
    a = modmul(s, kdepth=0, wa=256, wb=256)
    b = modmul(s, kdepth=0, wa=256, wb=513)
    print(f"    256 x 256 : {a['vars']:10,d} vars")
    print(f"    256 x 513 : {b['vars']:10,d} vars   ({b['vars'] / a['vars']:.2f}x)")
