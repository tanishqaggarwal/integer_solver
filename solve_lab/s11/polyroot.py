"""Root finding in GF(p)[x] for low-degree polynomials."""
import random


def norm(f, p):
    f = [c % p for c in f]
    while f and f[-1] == 0:
        f.pop()
    return f


def pmul(a, b, p):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
    return norm(r, p)


def pmod(a, b, p):
    a = norm(a[:], p)
    b = norm(b[:], p)
    if not b:
        raise ZeroDivisionError
    inv = pow(b[-1], -1, p)
    while len(a) >= len(b):
        c = a[-1] * inv % p
        sh = len(a) - len(b)
        for i, x in enumerate(b):
            a[i + sh] = (a[i + sh] - c * x) % p
        a = norm(a, p)
    return a


def pgcd(a, b, p):
    a = norm(a[:], p)
    b = norm(b[:], p)
    while b:
        a, b = b, pmod(a, b, p)
    if a:
        inv = pow(a[-1], -1, p)
        a = [c * inv % p for c in a]
    return a


def powmod(base, e, mod, p):
    r = [1]
    base = pmod(base, mod, p)
    while e:
        if e & 1:
            r = pmod(pmul(r, base, p), mod, p)
        base = pmod(pmul(base, base, p), mod, p)
        e >>= 1
    return r


def roots(f, p, rng=None):
    """all roots in GF(p) of f"""
    f = norm(f[:], p)
    if not f:
        raise ValueError("zero polynomial")
    if len(f) == 1:
        return []
    # g = gcd(x^p - x, f) = product of distinct linear factors
    xp = powmod([0, 1], p, f, p)
    g = pgcd(norm([(xp[0] if xp else 0)] + ([(xp[1] - 1) % p] + xp[2:] if len(xp) > 1 else [p - 1]), p)
             if False else norm([(-1 if i == 1 else 0) + (xp[i] if i < len(xp) else 0) for i in range(max(2, len(xp)))], p),
             f, p)
    out = []
    rng = rng or random.Random(0)
    stack = [g]
    while stack:
        h = norm(stack.pop(), p)
        if len(h) <= 1:
            continue
        if len(h) == 2:
            out.append((-h[0]) * pow(h[1], -1, p) % p)
            continue
        # equal-degree split
        for _ in range(200):
            a = rng.randrange(p)
            t = powmod([a, 1], (p - 1) // 2, h, p)
            t = norm([(t[0] - 1) % p] + t[1:], p) if t else norm([p - 1], p)
            d = pgcd(t, h, p)
            if 0 < len(d) - 1 < len(h) - 1:
                stack.append(d)
                stack.append(pdiv(h, d, p))
                break
        else:
            raise RuntimeError("split failed")
    return sorted(set(out))


def pdiv(a, b, p):
    a = norm(a[:], p)
    b = norm(b[:], p)
    q = [0] * max(1, len(a) - len(b) + 1)
    inv = pow(b[-1], -1, p)
    while len(a) >= len(b) and a:
        c = a[-1] * inv % p
        sh = len(a) - len(b)
        q[sh] = c
        for i, x in enumerate(b):
            a[i + sh] = (a[i + sh] - c * x) % p
        a = norm(a, p)
    return norm(q, p)


if __name__ == '__main__':
    P = 2**256 - 2**32 - 977
    rng = random.Random(1)
    for trial in range(5):
        rs = [rng.randrange(P) for _ in range(3)]
        f = [1]
        for r in rs:
            f = pmul(f, [(-r) % P, 1], P)
        got = roots(f, P, rng)
        print(trial, sorted(set(rs)) == got)
