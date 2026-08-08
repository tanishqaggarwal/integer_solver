"""ecsmall.py -- toy curves y^2 = x^3 + B over F_p, same shape as the real instance."""
def curve(p, B):
    def add(P, Q):
        if P is None: return Q
        if Q is None: return P
        x1, y1 = P; x2, y2 = Q
        if x1 == x2 and (y1 + y2) % p == 0: return None
        lam = (3*x1*x1 % p * pow(2*y1, -1, p) if P == Q else (y2-y1) * pow(x2-x1, -1, p)) % p
        x3 = (lam*lam - x1 - x2) % p
        return (x3, (lam*(x1-x3) - y1) % p)
    def mul(k, P):
        R = None
        while k:
            if k & 1: R = add(R, P)
            P = add(P, P); k >>= 1
        return R
    return add, mul

def find(p, B):
    add, mul = curve(p, B)
    pts = [(x, y) for x in range(p) for y in range(p) if (y*y - x*x*x - B) % p == 0]
    best = None
    for P in pts:
        o, Q = 1, P
        while Q is not None:
            Q = add(Q, P); o += 1
        if best is None or o > best[1]: best = (P, o)
    return best
