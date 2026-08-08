#!/usr/bin/env python3
"""gen.py -- synthetic ECDLP instances with a PLANTED, KNOWN key.

Same shape as the real instance (y^2 = x^3 + B over F_p, a doubling-chain comb,
a target T = k*G) but on a curve of chosen bit-size, with k generated here and
written alongside the instance.  No live target: every instance is solvable and
checkable end to end.

    inst = make(bits=32)          # 32-bit prime-order curve, random planted k
    inst.k                        # the planted scalar (the "answer")
    inst.verify()                 # k*G == T, on the nose

A curve is chosen with (near-)prime group order so the scheme matches the real
one (prime order => no Pohlig-Hellman shortcut, the honest setting).
"""
import random
from dataclasses import dataclass


def _is_prime(m):
    if m < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if m % q == 0: return m == q
    d, r = m - 1, 0
    while d % 2 == 0: d //= 2; r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, m); 
        if x in (1, m - 1): continue
        for _ in range(r - 1):
            x = x * x % m
            if x == m - 1: break
        else: return False
    return True


def _next_prime(m):
    m |= 1
    while not _is_prime(m): m += 2
    return m


class Curve:
    def __init__(self, p, B):
        self.p, self.B = p, B
    def on(self, P):
        if P is None: return True
        x, y = P; return (y*y - x*x*x - self.B) % self.p == 0
    def add(self, P, Q):
        p = self.p
        if P is None: return Q
        if Q is None: return P
        x1, y1 = P; x2, y2 = Q
        if x1 == x2 and (y1 + y2) % p == 0: return None
        lam = (3*x1*x1 % p * pow(2*y1, -1, p) if P == Q
               else (y2 - y1) * pow(x2 - x1, -1, p)) % p
        x3 = (lam*lam - x1 - x2) % p
        return (x3, (lam*(x1 - x3) - y1) % p)
    def mul(self, k, P):
        if k < 0: return self.mul(-k, (P[0], (-P[1]) % self.p)) if P else None
        R = None
        while k:
            if k & 1: R = self.add(R, P)
            P = self.add(P, P); k >>= 1
        return R
    def order_mult(self, P):
        """a multiple M of ord(P) in the Hasse interval, via BSGS. #E if that is prime."""
        import math
        p = self.p
        w = 1 + int((4 * math.isqrt(p)) ** 0.5)          # giant step
        # baby steps: j*P for j in [0, w]
        baby = {}
        R = None
        for j in range(w + 1):
            key = None if R is None else R[0]
            baby.setdefault(key, j); R = self.add(R, P)
        wP = self.mul(w, P)
        # Q0 = (p+1)P ; search (p+1 + a*w + b) P = O, |a*w+b| <= 2 sqrt p
        Q0 = self.mul((p + 1) % (0 or 10**30), P) if False else self.mul(p + 1, P)
        lim = 2 * math.isqrt(p) + 1
        A = w
        a = 0
        # walk Q = Q0 + a*wP for a = 0, -1, +1, ... looking for -bP among baby
        from itertools import count
        seen = {}
        cur = Q0
        step = wP
        cands = []
        amax = lim // w + 2
        for a in range(-amax, amax + 1):
            Qa = self.add(Q0, self.mul(a, wP))
            key = None if Qa is None else Qa[0]
            if key in baby:
                for b in (baby[key], -baby[key]):
                    # Qa == +-bP  => Q0 + a w P = -+ b P (x matches, check exactly)
                    if self.mul(p + 1 + a*w - b, P) is None:
                        cands.append(p + 1 + a*w - b)
                    if self.mul(p + 1 + a*w + b, P) is None:
                        cands.append(p + 1 + a*w + b)
        cands = [m for m in cands if m > 0]
        return min(cands) if cands else None


@dataclass
class Instance:
    curve: Curve
    G: tuple
    n: int          # order of G (prime)
    k: int          # PLANTED scalar
    T: tuple        # k*G
    bits: int
    def pts(self, m=None):
        """the doubling chain P_i = 2^i G, i = 0..bits-1 (or m-1)."""
        m = m or self.bits
        out, cur = [], self.G
        for _ in range(m): out.append(cur); cur = self.curve.add(cur, cur)
        return out
    def verify(self):
        c = self.curve
        assert c.on(self.G) and c.on(self.T)
        assert c.mul(self.n, self.G) is None
        assert c.mul(self.k, self.G) == self.T
        return True


def make(bits, seed=0, want_prime_order=True):
    """A prime-field curve whose base point has prime order of ~`bits` bits."""
    rnd = random.Random((seed << 20) ^ bits)
    for _try in range(2000):
        p = _next_prime(rnd.randrange(1 << (bits - 1), 1 << bits))
        B = rnd.randrange(1, p)
        c = Curve(p, B)
        # a random point
        for _ in range(40):
            x = rnd.randrange(p)
            rhs = (x*x*x + B) % p
            if pow(rhs, (p - 1) // 2, p) != 1: continue
            y = pow(rhs, (p + 1) // 4, p) if p % 4 == 3 else _tonelli(rhs, p)
            if y is None: continue
            P = (x, y)
            M = c.order_mult(P)
            if M is None: continue
            # accept only prime group order: then ord(P) = M and the group is cyclic prime
            if not _is_prime(M): continue
            if c.mul(M, P) is not None: continue
            o = M
            if o < (1 << (bits - 2)): continue
            k = rnd.randrange(2, o)
            inst = Instance(c, P, o, k, c.mul(k, P), bits)
            inst.verify()
            return inst
    raise RuntimeError(f"no {bits}-bit prime-order curve found")


def _tonelli(a, p):
    if pow(a, (p - 1) // 2, p) != 1: return None
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    if s == 1: return pow(a, (p + 1) // 4, p)
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, cc, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1: tt = tt * tt % p; i += 1
        b = pow(cc, 1 << (m - i - 1), p)
        m, cc, t, r = i, b*b % p, t*b*b % p, r*b % p
    return r


if __name__ == '__main__':
    for bits in (8, 12, 16, 24, 32, 40):
        inst = make(bits)
        print(f"bits={bits:3d}  p={inst.curve.p:<14} B={inst.curve.B:<14} "
              f"n={inst.n:<14} (prime={_is_prime(inst.n)})  k={inst.k}  verified={inst.verify()}")
