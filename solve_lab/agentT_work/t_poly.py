"""Modular univariate root-finding for AUDIT T31/T32 -- needed because the two-wire root search
cannot brute-force q^2 when c has a large prime factor.  Degrees here are <= 5.
Everything is checked against brute force in the self-test at the bottom."""
import random
from math import gcd as _g

def strip(f):
    while f and f[-1] == 0:
        f.pop()
    return f

def pmul(a, b, q):
    if not a or not b:
        return []
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i+j] = (r[i+j] + x*y) % q
    return strip(r)

def psub(a, b, q):
    n = max(len(a), len(b)); r = [0]*n
    for i in range(n):
        r[i] = ((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % q
    return strip(r)

def pmod(a, b, q):
    a = a[:]; db = len(b)-1
    if db < 0:
        raise ZeroDivisionError
    ib = pow(b[-1], -1, q)
    while len(a)-1 >= db and a:
        d = len(a)-1-db
        c = a[-1]*ib % q
        for i in range(db+1):
            a[d+i] = (a[d+i]-c*b[i]) % q
        strip(a)
    return a

def pdiv(a, b, q):
    a = a[:]; db = len(b)-1; ib = pow(b[-1], -1, q); out = [0]*max(0, len(a)-db)
    while len(a)-1 >= db and a:
        d = len(a)-1-db
        c = a[-1]*ib % q
        out[d] = c
        for i in range(db+1):
            a[d+i] = (a[d+i]-c*b[i]) % q
        strip(a)
    return strip(out)

def pgcd(a, b, q):
    a = strip(a[:]); b = strip(b[:])
    while b:
        a, b = b, pmod(a, b, q)
    if a:
        ia = pow(a[-1], -1, q)
        a = [x*ia % q for x in a]
    return a

def ppow(base, e, mod, q):
    r = [1]; base = pmod(base[:], mod, q)
    while e:
        if e & 1:
            r = pmod(pmul(r, base, q), mod, q)
        base = pmod(pmul(base, base, q), mod, q)
        e >>= 1
    return r

def _split(g, q):
    g = strip(g[:])
    d = len(g)-1
    if d <= 0:
        return []
    ig = pow(g[-1], -1, q); g = [x*ig % q for x in g]
    if d == 1:
        return [(-g[0]) % q]
    while True:
        a = random.randrange(q)
        h = psub(ppow([a, 1], (q-1)//2, g, q), [1], q)
        d1 = pgcd(h, g, q)
        if 0 < len(d1)-1 < d:
            return _split(d1, q) + _split(pdiv(g, d1, q), q)

def roots_prime(f, q, brute=3000):
    """all roots mod prime q of monomial-coefficient poly f (low->high).  [] if none.
    returns the string 'ALL' if f is identically zero mod q."""
    f = strip([x % q for x in f])
    if not f:
        return 'ALL'
    if q <= brute or len(f)-1 <= 1:
        if len(f)-1 == 1:
            return [(-f[0])*pow(f[1], -1, q) % q]
        out = []
        for t in range(q):
            v = 0
            for c in reversed(f):
                v = (v*t + c) % q
            if v == 0:
                out.append(t)
        return out
    if q == 2:
        return [t for t in (0, 1) if sum(c*(t**i) for i, c in enumerate(f)) % 2 == 0]
    xq = ppow([0, 1], q, f, q)
    g = pgcd(psub(xq, [0, 1], q), f, q)
    return sorted(_split(g, q))

def roots_pp(f, m, q, e, brute=3000):
    """roots mod q^e, by Hensel/enumeration lifting from roots mod q."""
    if e == 1:
        return roots_prime(f, q, brute)
    if m <= brute:
        out = []
        for t in range(m):
            v = 0
            for c in reversed(f):
                v = (v*t + c) % m
            if v == 0:
                out.append(t)
        return out
    r = roots_prime(f, q, brute)
    if r == 'ALL':
        return 'ALL'
    df = [i*c for i, c in enumerate(f)][1:]
    cur = r; mod = q
    for _ in range(e-1):
        nxt = []
        for t in cur:
            for k in range(q):
                tt = t + mod*k
                v = 0
                for c in reversed(f):
                    v = (v*tt + c) % (mod*q)
                if v == 0:
                    nxt.append(tt)
        cur = sorted(set(nxt)); mod *= q
    return cur

def newton_to_mono(cf, m):
    """Newton-basis coeffs (sum_k cf[k]*C(t,k)) -> monomial coeffs mod m.  Requires k! invertible
    mod m for every k with cf[k] != 0; returns None if it is not."""
    out = [0]*len(cf)
    ff = [1]                       # falling factorial prod_{j<k}(t-j), integer coeffs
    fact = 1
    for k in range(len(cf)):
        if k:
            ff = [(ff[i-1] if i else 0) - (k-1)*(ff[i] if i < len(ff) else 0)
                  for i in range(len(ff)+1)]
            fact *= k
        if cf[k] % m:
            if _g(fact, m) != 1:
                return None
            inv = pow(fact, -1, m)
            c = cf[k] % m * inv % m
            for i, x in enumerate(ff):
                out[i] = (out[i] + c*x) % m
    return strip(out)

if __name__ == '__main__':
    random.seed(3)
    ok = bad = 0
    for _ in range(300):
        q = random.choice([3, 7, 19, 83, 463, 4001, 15013, 104729])
        d = random.randrange(1, 6)
        f = [random.randrange(q) for _ in range(d+1)]
        if f[-1] == 0:
            f[-1] = 1
        got = roots_prime(f, q)
        want = [t for t in range(q) if sum(c*pow(t, i, q) for i, c in enumerate(f)) % q == 0] \
            if q <= 200000 else None
        if want is not None and sorted(got) != want:
            bad += 1; print('MISMATCH', q, f, sorted(got), want)
        else:
            ok += 1
    print('roots_prime self-test: %d ok, %d mismatched' % (ok, bad))
    # newton_to_mono self-test
    bad2 = 0
    for _ in range(200):
        m = random.choice([15013, 104729, 4001])   # prime > deg, so k! is invertible
        cf = [random.randrange(-10**6, 10**6) for _ in range(6)]
        mo = newton_to_mono(cf, m)
        for t in range(5):
            a = 0; b = 1
            for k in range(6):
                if k:
                    b = b*(t-k+1)//k
                a += cf[k]*b
            v = sum(c*pow(t, i, m) for i, c in enumerate(mo)) % m
            if a % m != v:
                bad2 += 1
    print('newton_to_mono self-test: %d mismatches' % bad2)
