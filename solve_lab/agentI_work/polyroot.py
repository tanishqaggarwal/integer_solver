"""Root finding for small polynomials over F_p."""
import random
P = 2**256 - 2**32 - 977


def trim(f):
    while f and f[-1] == 0:
        f.pop()
    return f


def pmul(a, b):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % P
    return trim(r)


def padd(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i in range(n):
        r[i] = ((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % P
    return trim(r)


def psub(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i in range(n):
        r[i] = ((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % P
    return trim(r)


def pmod(a, m):
    a = a[:]
    dm = len(m) - 1
    inv = pow(m[-1], -1, P)
    while len(a) - 1 >= dm and a:
        d = len(a) - 1 - dm
        c = a[-1] * inv % P
        for i in range(len(m)):
            a[d + i] = (a[d + i] - c * m[i]) % P
        trim(a)
    return a


def pgcd(a, b):
    a = trim(a[:]); b = trim(b[:])
    while b:
        a, b = b, pmod(a, b)
    if a:
        inv = pow(a[-1], -1, P)
        a = [x * inv % P for x in a]
    return a


def powmod(base, e, m):
    r = [1]; b = pmod(base, m)
    while e:
        if e & 1:
            r = pmod(pmul(r, b), m)
        b = pmod(pmul(b, b), m)
        e >>= 1
    return r


def roots(f):
    f = trim(f[:])
    if not f:
        return None    # zero polynomial
    if len(f) == 1:
        return []
    xp = powmod([0, 1], P, f)
    g = pgcd(psub(xp, [0, 1]), f)
    out = []
    _split(g, out)
    return sorted(set(out))


def _split(g, out):
    g = trim(g[:])
    if len(g) <= 1:
        return
    if len(g) == 2:
        out.append((-g[0]) * pow(g[1], -1, P) % P)
        return
    rng = random.Random(12345)
    while True:
        a = rng.randrange(P)
        h = powmod([a, 1], (P - 1) // 2, g)
        d = pgcd(psub(h, [1]), g)
        if 1 <= len(d) - 1 < len(g) - 1:
            _split(d, out)
            _split(pdiv(g, d), out)
            return


def pdiv(a, m):
    a = a[:]
    q = [0] * max(1, len(a) - len(m) + 1)
    dm = len(m) - 1
    inv = pow(m[-1], -1, P)
    while len(a) - 1 >= dm and a:
        d = len(a) - 1 - dm
        c = a[-1] * inv % P
        q[d] = c
        for i in range(len(m)):
            a[d + i] = (a[d + i] - c * m[i]) % P
        trim(a)
    return trim(q)
