"""Agent AD -- small analogue, step 1: find small curves with prime order.

Short Weierstrass y^2 = x^3 + a*x + b over F_p, p prime.
We want  N = #E(F_p)  PRIME and of exactly `nbits` bits, because the real
instance has prime order N (secp256k1) and n = 256 leaf selectors = bitlen(N).

At least one j=0 curve (a=0, CM by sqrt(-3)) per size, because the real
instance is j=0.  Also non-j=0 curves so we can see whether j=0 matters.

Output: ad_curves.json
"""
import json, os, sys
from sympy import isprime

HERE = os.path.dirname(os.path.abspath(__file__))


def qr_table(p):
    t = bytearray(p)
    for x in range((p + 1) // 2 + 1):
        t[(x * x) % p] = 1
    return t


def order_of_curve(p, a, b, qr):
    # #E = p + 1 + sum_x chi(f(x)) ; chi(0)=0, chi(qr)=1, chi(nqr)=-1
    tot = p + 1
    for x in range(p):
        v = (x * x % p * x + a * x + b) % p
        if v == 0:
            continue
        tot += 1 if qr[v] else -1
    return tot


def point_add(P1, P2, p, a):
    if P1 is None:
        return P2
    if P2 is None:
        return P1
    x1, y1 = P1
    x2, y2 = P2
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1 + a) * pow(2 * y1 % p, -1, p) % p
    else:
        lam = (y2 - y1) * pow((x2 - x1) % p, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def scalar_mul(k, P, p, a):
    R = None
    Q = P
    while k:
        if k & 1:
            R = point_add(R, Q, p, a)
        Q = point_add(Q, Q, p, a)
        k >>= 1
    return R


def find_gen(p, a, b, N):
    for x in range(1, p):
        v = (x * x % p * x + a * x + b) % p
        y = pow(v, (p + 1) // 4, p) if p % 4 == 3 else tonelli(v, p)
        if y is None:
            continue
        if (y * y - v) % p:
            continue
        P = (x, y)
        if scalar_mul(N, P, p, a) is None and P is not None:
            return P
    return None


def tonelli(n, p):
    n %= p
    if n == 0:
        return 0
    if pow(n, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        bb = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, bb * bb % p, t * bb % p * bb % p, r * bb % p
    return r


def primes_near(nbits, count, cong3=None):
    out = []
    v = (1 << nbits) - 1
    while v > (1 << (nbits - 1)) and len(out) < count:
        if isprime(v) and (cong3 is None or v % 3 == cong3):
            out.append(v)
        v -= 2
    return out


def search(nbits, want_j0=2, want_gen=3, seed=1):
    """Return list of curve dicts with prime order of exactly nbits bits."""
    import random
    rng = random.Random(seed + nbits)
    found = []
    # j=0 curves need p = 1 mod 3 for the sextic-twist family to be rich
    for p in primes_near(nbits, 8, cong3=1):
        qr = qr_table(p)
        # ---- j = 0 : a = 0
        n0 = 0
        bs = list(range(1, min(p, 400)))
        rng.shuffle(bs)
        for b in bs:
            if n0 >= want_j0:
                break
            N = order_of_curve(p, 0, b, qr)
            if isprime(N) and N.bit_length() == nbits:
                found.append(dict(p=p, a=0, b=b, N=N, j0=True))
                n0 += 1
        # ---- generic  a != 0
        ng = 0
        tries = 0
        while ng < want_gen and tries < 400:
            tries += 1
            a = rng.randrange(1, p)
            b = rng.randrange(0, p)
            if (4 * a ** 3 + 27 * b * b) % p == 0:
                continue
            N = order_of_curve(p, a, b, qr)
            if isprime(N) and N.bit_length() == nbits:
                found.append(dict(p=p, a=a, b=b, N=N, j0=False))
                ng += 1
        if sum(1 for c in found if c['j0']) >= want_j0 and \
           sum(1 for c in found if not c['j0']) >= want_gen:
            break
    return found


def main():
    out = {}
    for nbits in (8, 12, 16, 20):
        cur = search(nbits)
        for c in cur:
            G = find_gen(c['p'], c['a'], c['b'], c['N'])
            assert G is not None
            assert scalar_mul(c['N'], G, c['p'], c['a']) is None
            c['G'] = list(G)
        out[str(nbits)] = cur
        print(nbits, 'bits:', len(cur), 'curves')
        for c in cur:
            print('   p=%d a=%d b=%d N=%d j0=%s G=%s' %
                  (c['p'], c['a'], c['b'], c['N'], c['j0'], c['G']))
        sys.stdout.flush()
    with open(os.path.join(HERE, 'ad_curves.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print('written')


if __name__ == '__main__':
    main()
