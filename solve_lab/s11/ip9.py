"""IP #9 -- WHAT divisibility fails?

The checkpoint's system  M x = rhs  is solvable mod every prime tried but not over Z.  The
obstruction is therefore an invariant-factor / denominator condition.  Solve over Q and report
the least common denominator D: that single integer IS the obstruction, and if it is small it
is attackable exactly the way the 8640431 condition was.
"""
import sys, os, json, time
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip8 import build
from ip7 import load_raw
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)


def rational_solve(M, rhs):
    """least-denominator rational solution of M x = rhs (or None if inconsistent over Q)"""
    m = len(M)
    n = len(M[0]) if m else 0
    A = [[Fraction(M[i][j]) for j in range(n)] + [Fraction(rhs[i])] for i in range(m)]
    piv = []
    r = 0
    for c in range(n):
        pr = None
        for i in range(r, m):
            if A[i][c] != 0:
                pr = i
                break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        f = A[r][c]
        A[r] = [x / f for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                g = A[i][c]
                A[i] = [A[i][k] - g * A[r][k] for k in range(n + 1)]
        piv.append(c)
        r += 1
        if r == m:
            break
    for i in range(r, m):
        if A[i][n] != 0 and all(A[i][j] == 0 for j in range(n)):
            return None, None
    x = [Fraction(0)] * n
    for i, c in enumerate(piv):
        x[c] = A[i][n]
    return x, piv


if __name__ == '__main__':
    LAB = os.path.join(HERE, '..')
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
    v = load_raw(src)
    print(f"=== {os.path.basename(src)}")
    t0 = time.time()
    v, FAIL, used, M, rhs, nf = build(v)
    print(f"  solving {len(M)} x {len(M[0])} over Q ...", flush=True)
    x, piv = rational_solve(M, rhs)
    if x is None:
        print("  INCONSISTENT even over Q -> a genuine rank obstruction, not divisibility")
        sys.exit()
    D = 1
    for t in x:
        D = D * t.denominator // __import__('math').gcd(D, t.denominator)
    print(f"  consistent over Q; least common denominator D has {len(str(D))} digits "
          f"({time.time()-t0:.0f}s)")
    print(f"  D = {str(D)[:120]}{'...' if len(str(D)) > 120 else ''}")
    # factor the small part of D
    d = D
    small = []
    for q in range(2, 200000):
        if q * q > d:
            break
        while d % q == 0:
            small.append(q)
            d //= q
    print(f"  small prime factors of D: {small[:20]}")
    print(f"  remaining cofactor: {len(str(d))} digits")
    nz = [i for i, t in enumerate(x) if t != 0]
    print(f"  rational solution support: {len(nz)} of {len(x)} variables")
    print(f"  variables used: {[used[i] for i in nz][:20]}")
