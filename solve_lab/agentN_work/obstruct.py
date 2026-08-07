"""Why is the integrally-zeroable row count pinned?  Localise the obstruction.

For a placement's region (M, b): the rows are rationally consistent for every subset
(measured), so the whole obstruction is integrality.  For each row subset test which prime
blocks it: compute the smallest d>0 with d*(-b_I) in the integer column lattice of M_I.
Every prime dividing d is an obstructing prime for I.
"""
import sys, json, itertools, time
from flint import fmpz_mat
import optN
from optN import make, build, WIT, POOL
import zsolve
from fractions import Fraction


def rational_solvable(M, b, rows, n):
    A = [[Fraction(M[i][j]) for j in range(n)] + [Fraction(-b[i])] for i in rows]
    r = 0
    for c in range(n):
        k = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if k is None:
            continue
        A[r], A[k] = A[k], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(n + 1)]
        r += 1
        if r == len(A):
            break
    for i in range(r, len(A)):
        if A[i][n] != 0 and all(A[i][j] == 0 for j in range(n)):
            return False
    return True


def denom(M, b, rows, n):
    """Smallest d>0 with d*(-b_I) in the Z-column-lattice of M_I; None if not even rational."""
    key = sorted(rows)
    gens = [[M[i][j] for i in key] for j in range(n)]
    H = zsolve._hnf_rows(gens)
    c = [-b[i] for i in key]
    m = len(c)
    piv = []
    for r in H:
        j = next((k for k in range(m) if r[k]), None)
        if j is not None:
            piv.append((j, r))
    piv.sort(key=lambda z: z[0])
    # solve over Q along the HNF, collecting denominators
    from fractions import Fraction as F
    cc = [F(x) for x in c]
    d = 1
    for j, r in piv:
        if cc[j] != 0:
            f = cc[j] / r[j]
            d = d * f.denominator // __import__('math').gcd(d, f.denominator)
            for k in range(j, m):
                cc[k] -= f * r[k]
    if any(x != 0 for x in cc):
        return None            # not in the rational span at all
    return d


def report(tag, D, kind='detach'):
    if kind == 'detach':
        st = make(list(D))
    d = build(st)
    M, b, n, Rl = d['M'], d['b'], d['n'], d['R']
    nr = len(Rl)
    print('\n=== %s  |R|=%d knobs=%d ===' % (tag, nr, n), flush=True)
    print('rationally solvable, ALL %d rows at once: %s' % (nr, rational_solvable(M, b, range(nr), n)),
          flush=True)
    opt, rows, exh, tests = zsolve.max_zero_rows(M, b, n, nr)
    print('integer optimum = %d  (rows %s)' % (opt, rows), flush=True)
    # obstruction primes for every subset of size opt+1
    from collections import Counter
    cnt = Counter()
    nsub = 0
    unsolv = 0
    for I in itertools.combinations(range(nr), opt + 1):
        nsub += 1
        dd = denom(M, b, I, n)
        if dd is None:
            cnt['NOT-RATIONAL'] += 1
            continue
        if dd == 1:
            continue
        unsolv += 1
        x = dd
        f = []
        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
                  73, 79, 83, 89, 97, 101, 103, 107, 109, 113):
            while x % p == 0:
                f.append(p)
                x //= p
        if x != 1:
            f.append(x)
        cnt[tuple(sorted(set(f)))] += 1
    print('subsets of size %d: %d, all blocked; obstruction-denominator prime signatures:' %
          (opt + 1, nsub), flush=True)
    for k, v in cnt.most_common(12):
        kk = [str(z)[:24] + ('..' if len(str(z)) > 24 else '') for z in k] if isinstance(k, tuple) else k
        print('   %-6d x  %s' % (v, kk), flush=True)
    return M, b, n, Rl


if __name__ == '__main__':
    P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
    M, b, n, Rl = report('WITNESS %s' % WIT, WIT)
    print('\np (the 256-bit constant in the lab notes) divides denominators?')
    import itertools as it
    hits = 0
    tot = 0
    for I in it.combinations(range(len(Rl)), 6):
        dd = denom(M, b, I, n)
        if dd and dd > 1:
            tot += 1
            if dd % P == 0:
                hits += 1
    print('  size-6 subsets blocked: %d, of which denominator divisible by p: %d' % (tot, hits))
    report('SINGLETON [17499]', [17499])
