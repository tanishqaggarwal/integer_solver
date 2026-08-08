#!/usr/bin/env python3
"""instance.py -- the verified decision core, as plain arithmetic. Import this.

    from instance import p, n, A, B, G, T, PTS, add, mul, order_ok

Everything here was derived in reduce.py / structure.py directly from the constants
in EQUATIONS.txt and re-checked on import.  Nothing is assumed.

    curve   y^2 = x^3 + B   over F_p,   p = 2^256 - 2^32 - 977   (A == 0)
    PTS[i]  = 2^i * G  for i = 0..255   (verified doubling chain)
    problem find b_0..b_255 in {0,1} with  sum b_i PTS[i] == T
            i.e.  k*G == T  with k = sum b_i 2^i
    n       = group order, PRIME
"""
import json, os

_H = os.path.dirname(os.path.abspath(__file__))
_C = json.load(open(os.path.join(_H, 'core.json')))
_CH = json.load(open(os.path.join(_H, 'chain.json')))

p = int(_C['p'])
A = int(_C['A'])                     # == 0
B = int(_C['B'])
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141   # prime


def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0: return None
    lam = (3*x1*x1 % p * pow(2*y1, -1, p) if P == Q else (y2-y1) * pow(x2-x1, -1, p)) % p
    x3 = (lam*lam - x1 - x2) % p
    return (x3, (lam*(x1-x3) - y1) % p)


def neg(P): return None if P is None else (P[0], (-P[1]) % p)
def sub(P, Q): return add(P, neg(Q))


def mul(k, P):
    if k < 0: return mul(-k, neg(P))
    R = None
    k %= n
    while k:
        if k & 1: R = add(R, P)
        P = add(P, P); k >>= 1
    return R


# the doubling chain, in order:  PTS[i] = 2^i * G
PTS = [tuple(int(v) for v in _C['points_dep'][str(b)]) for b in _CH['chain_bit_vars']]
G = PTS[0]
T = tuple(int(v) for v in _C['target_dep'])

# --- the F_p-isomorphism to  y^2 = x^3 + 7  (secp256k1 short form) ---
# (x, y)  ->  (x/u^2, y/u^3)   with  u^6 = B/7
def _sixth_root():
    from sympy.ntheory.residue_ntheory import nthroot_mod
    return nthroot_mod(B * pow(7, -1, p) % p, 6, p)


U = None
def to_secp(P):
    """map a point of this curve onto y^2 = x^3 + 7."""
    global U
    if U is None: U = _sixth_root()
    if P is None: return None
    return (P[0] * pow(U, -2, p) % p, P[1] * pow(U, -3, p) % p)


def from_secp(P):
    global U
    if U is None: U = _sixth_root()
    if P is None: return None
    return (P[0] * pow(U, 2, p) % p, P[1] * pow(U, 3, p) % p)


SECP_G = (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
          0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)


def selfcheck():
    assert A == 0
    assert all((q[1]*q[1] - q[0]**3 - B) % p == 0 for q in PTS), "point off curve"
    assert (T[1]*T[1] - T[0]**3 - B) % p == 0, "target off curve"
    assert mul(n, G) is None and mul(n, T) is None, "order wrong"
    for i in range(1, 256):
        assert PTS[i] == add(PTS[i-1], PTS[i-1]), f"chain break at {i}"
    return True


if __name__ == '__main__':
    selfcheck()
    print("selfcheck OK")
    print("p =", p)
    print("B =", B)
    print("n =", n)
    print("G =", G)
    print("T =", T)
    print("G on secp256k1 coords =", to_secp(G))
    print("T on secp256k1 coords =", to_secp(T))
    print("is G the standard secp256k1 generator?", to_secp(G) == SECP_G)
