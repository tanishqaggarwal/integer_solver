#!/usr/bin/env python3
"""agent V -- multivariate exact fit + root finding mod q^e, for ARBITRARY arity k.

Written from scratch rather than borrowed, so it is a second implementation of the same idea
rather than a second run of the same code.  Every consumer verifies its answers by DIRECT
RECOMPUTATION, so a defect here can only cause a MISSED solution, never a false one -- but the
fit is additionally validated against recomputation at random points before it is used.

Representation: a k-dimensional Newton forward-difference array `cf`, indexed cf[i1]...[ik], with
    R(t1..tk)/p  =  sum_{i}  cf[i] * prod_j C(tj, ij)         (exact over Z)
Newton form is used because C(t,i) has integer values at integer t and the forward-difference
fit is exact with no division.
"""
import itertools
from math import gcd


# ------------------------------------------------------------------ Newton fit, arity k
def nd_fit(probe_fn, k, D):
    """probe_fn(tuple of k ints) -> integer value of R/p.  Returns the k-dim Newton array."""
    shape = (D+1,)*k
    grid = {}
    for idx in itertools.product(range(D+1), repeat=k):
        grid[idx] = probe_fn(idx)
    cur = grid
    for axis in range(k):
        nxt = {}
        for idx in itertools.product(*[range(D+1) for _ in range(k)]):
            # forward differences along `axis`
            pass
        # do it explicitly: for each fixed value of the other axes, difference along `axis`
        others = [range(D+1) for _ in range(k)]
        others[axis] = [None]
        for pre in itertools.product(*others):
            row = []
            for t in range(D+1):
                idx = list(pre); idx[axis] = t
                row.append(cur[tuple(idx)])
            d = [row[:]]
            for s in range(D):
                d.append([d[s][j+1]-d[s][j] for j in range(len(d[s])-1)])
            for s in range(D+1):
                idx = list(pre); idx[axis] = s
                nxt[tuple(idx)] = d[s][0]
        cur = nxt
    return cur                          # dict: multi-index -> Newton coefficient


def nd_eval(cf, ts, m=None):
    """evaluate the Newton array at integer point ts; exactly if m is None, else mod m."""
    k = len(ts)
    D = max(i[0] for i in cf) if cf else 0
    for idx in cf:
        D = max(D, max(idx))
    B = []
    for t in ts:
        b = [1]
        for i in range(1, D+1):
            b.append(b[-1]*(t-i+1)//i)
        B.append(b)
    tot = 0
    for idx, c in cf.items():
        if not c:
            continue
        v = c
        for j, i in enumerate(idx):
            v *= B[j][i]
        tot += v
    return tot if m is None else tot % m


def nd_degrees(cf):
    """observed top degree per variable (a nonzero Newton coefficient at that index)."""
    if not cf:
        return []
    k = len(next(iter(cf)))
    return [max((idx[j] for idx, c in cf.items() if c), default=0) for j in range(k)]


# ------------------------------------------------------------------ univariate, Newton -> monomial
def newton_to_mono(u, m):
    """convert Newton coefficients u (sum u_i C(t,i)) to monomial coefficients mod m.
    C(t,i) = t(t-1)...(t-i+1)/i!  -- needs i! invertible mod m, so only valid when m > deg and
    gcd(i!, m) == 1.  Returns None when that fails."""
    D = len(u)-1
    out = [0]*(D+1)
    for i in range(D+1):
        if not u[i]:
            continue
        f = 1
        for j in range(1, i+1):
            f *= j
        if gcd(f, m) != 1:
            return None
        # falling factorial t(t-1)...(t-i+1) as monomial coefficients
        poly = [1]
        for j in range(i):
            nxt = [0]*(len(poly)+1)
            for d, c in enumerate(poly):
                nxt[d+1] += c
                nxt[d] -= c*j
            poly = nxt
        inv = pow(f, -1, m)
        for d, c in enumerate(poly):
            out[d] = (out[d] + u[i]*c*inv) % m
    return out


def mono_eval(mo, t, m):
    v = 0
    for c in reversed(mo):
        v = (v*t + c) % m
    return v


def roots_prime(mo, q):
    """all roots mod prime q of a monomial-coefficient poly.  'ALL' if identically zero mod q."""
    mo = [c % q for c in mo]
    while mo and mo[-1] == 0:
        mo.pop()
    if not mo:
        return 'ALL'
    d = len(mo)-1
    if d == 0:
        return []
    if d == 1:
        a, b = mo[1], mo[0]
        g = gcd(a, q)
        if g == 1:
            return [(-b)*pow(a, -1, q) % q]
        return [t for t in range(q) if (a*t+b) % q == 0]
    if q <= 4096 or d*d*16 > q:
        return [t for t in range(q) if mono_eval(mo, t, q) == 0]
    # large prime, small degree: gcd(x^q - x, f) then split.  Degree <= 4 here, so simply
    # test the roots of the gcd by trial over its own (tiny) factorisation is not needed --
    # use Cantor-Zassenhaus-lite: find gcd(f, x^q-x) then equal-degree split.
    return _roots_large_prime(mo, q)


def _pmul(a, b, f, q):
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if not x:
            continue
        for j, y in enumerate(b):
            r[i+j] = (r[i+j] + x*y) % q
    return _pmod(r, f, q)


def _pmod(a, f, q):
    a = a[:]
    df = len(f)-1
    inv = pow(f[-1], -1, q)
    while len(a)-1 >= df and any(a):
        while a and a[-1] == 0:
            a.pop()
        if len(a)-1 < df:
            break
        c = a[-1]*inv % q
        sh = len(a)-1-df
        for i, x in enumerate(f):
            a[sh+i] = (a[sh+i] - c*x) % q
        while a and a[-1] == 0:
            a.pop()
    return a or [0]


def _pgcd(a, b, q):
    a = [x % q for x in a]; b = [x % q for x in b]
    while b and any(b):
        a, b = b, _pmod(a, b, q)
        while b and b[-1] == 0:
            b.pop()
    while a and a[-1] == 0:
        a.pop()
    return a or [0]


def _ppow(base, ex, f, q):
    r = [1]
    while ex:
        if ex & 1:
            r = _pmul(r, base, f, q)
        base = _pmul(base, base, f, q)
        ex >>= 1
    return r


def _pdiv(a, b, q):
    a = a[:]; db = len(b)-1; inv = pow(b[-1], -1, q)
    out2 = [0]*max(1, len(a)-db)
    while any(a):
        while a and a[-1] == 0:
            a.pop()
        if len(a)-1 < db:
            break
        c = a[-1]*inv % q
        sh = len(a)-1-db
        out2[sh] = c
        for i, x in enumerate(b):
            a[sh+i] = (a[sh+i] - c*x) % q
        while a and a[-1] == 0:
            a.pop()
    return out2


def _roots_large_prime(mo, q):
    import random
    f = [c % q for c in mo]
    while f and f[-1] == 0:
        f.pop()
    # g = gcd(f, x^q - x): the product of (x-r) over all roots r
    xq = _ppow([0, 1], q, f, q)
    xq = xq[:] + [0]*max(0, 2-len(xq))
    sub = xq[:]
    while len(sub) < 2:
        sub.append(0)
    sub[1] = (sub[1] - 1) % q
    g = _pgcd(f, sub, q)
    if len(g)-1 <= 0:
        return []
    out = []
    rnd = random.Random(1234567)

    def split(h):
        d = len(h)-1
        if d == 0:
            return
        if d == 1:
            out.append((-h[0])*pow(h[1], -1, q) % q)
            return
        while True:
            a = rnd.randrange(q)
            t = _ppow([a, 1], (q-1)//2, h, q)
            t = t[:] + [0]*max(0, 1-len(t))
            t[0] = (t[0]-1) % q
            g2 = _pgcd(h, t, q)
            if 0 < len(g2)-1 < d:
                split(g2)
                split(_pdiv(h, g2, q))
                return

    split(g)
    return sorted(set(out))


def roots_pp(u_newton, q, e):
    """roots mod q^e of the Newton-form univariate poly u.  Returns a list, or 'ALL'.
    Hensel-lifts the mod-q roots; falls back to enumeration when q^e is small."""
    m = q**e
    # MEASURED: this is the inner loop of the joint solve -- it runs q^(e(k-1)) times, so an
    # O(m) enumeration here costs O(m^k) overall.  At q=39703 that is 1.6e9 evaluations (~35 min
    # observed) versus ~8 s through the monomial route.  Enumerate only when m is genuinely tiny.
    if m <= 1024:
        return [t for t in range(m) if _neval(u_newton, t, m) == 0]
    if e == 1:
        mo = newton_to_mono(u_newton, q)
        if mo is None:
            return [t for t in range(q) if _neval(u_newton, t, q) == 0]
        r = roots_prime(mo, q)
        return list(range(q)) if r == 'ALL' else sorted(r)
    mo = newton_to_mono(u_newton, m)
    if mo is None:
        return [t for t in range(m) if _neval(u_newton, t, m) == 0]
    base = roots_prime([c % q for c in mo], q)
    if base == 'ALL':
        base = list(range(q))
    cur = base
    for lev in range(2, e+1):
        mm = q**lev
        nxt = []
        for r in cur:
            for j in range(q):
                t = r + j*q**(lev-1)
                if mono_eval(mo, t, mm) == 0:
                    nxt.append(t)
        cur = nxt
        if not cur:
            return []
    return sorted(set(cur))


def _neval(u, t, m):
    tot = 0; b = 1
    for i, c in enumerate(u):
        if i:
            b = b*(t-i+1)//i
        tot += c*b
    return tot % m


def mono_eval_exact(mo, t):
    v = 0
    for c in reversed(mo):
        v = v*t + c
    return v
