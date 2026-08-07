#!/usr/bin/env python3
"""Option 2 measured: bit-blast only the SELECTOR layer, treat stage values as an
uninterpreted binary operation with the law given only as axioms (assoc + comm).
Question: does the solver then search configurations usefully?"""
import time
from z3 import *

def run(n, extra_axioms=False, tmo=60):
    S = DeclareSort('S')
    f = Function('f', S, S, S)
    Ls = [Const('L%d' % i, S) for i in range(n)]
    T = Const('T', S)
    x, y, z = Consts('x y z', S)
    s = Solver(); s.set('timeout', tmo * 1000)
    s.add(ForAll([x, y], f(x, y) == f(y, x)))
    s.add(ForAll([x, y, z], f(f(x, y), z) == f(x, f(y, z))))
    if extra_axioms:
        s.add(Distinct(*(Ls + [T])))
    b = [Bool('b%d' % i) for i in range(n)]
    acc = Ls[0]
    for i in range(1, n):
        nxt = Const('a%d' % i, S)
        s.add(nxt == If(b[i], f(acc, Ls[i]), acc))
        acc = nxt
    s.add(acc == T)
    t = time.time()
    r = s.check()
    return r, round(time.time() - t, 2), (s.model() if r == sat else None)

if __name__ == '__main__':
    for n in (8, 16, 32, 64, 128, 256):
        for ax in (False, True):
            r, t, m = run(n, ax)
            print('n=%3d distinct-axiom=%-5s -> %-7s %6.2fs' % (n, ax, r, t), flush=True)
