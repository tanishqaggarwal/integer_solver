#!/usr/bin/env python3
"""The gap obstruction priced in tradeoff.py is an ARTEFACT: it assumed the intervening
selectors are arbitrary.  Zeroing them makes the mux the identity, acc does not move, and the
elimination stays low degree (collapse.py confirms solutions exist at every gap on siblings).

So solve it at the REAL 256-bit prime.  With base accumulator A and relaxed leaves i<j:
   acc = A + t1*(chord(A,L_i) - A)                      -- degree 1 in t1
   need T on the line through acc and chord(acc,L_j)    -- one poly in t1, degree ~6
Root-find that, then t2 follows linearly."""
import sys, json
sys.path.insert(0, '.')
import model
from model import P
lad = [(int(a), int(b)) for _, a, b in json.load(open('ladder.json'))['ladder']]
T = model.to_short(model.TARGET)

def chord(A, B):
    ax, ay = A; bx, by = B
    if (ax - bx) % P == 0: return None
    l = (by - ay) * pow(bx - ax, P - 2, P) % P
    sx = (l * l - ax - bx) % P
    return (sx, (l * (ax - sx) - ay) % P)

def resid(A, i, j, t1):
    Si = chord(A, lad[i])
    if Si is None: return None
    acc = ((A[0] + t1 * (Si[0] - A[0])) % P, (A[1] + t1 * (Si[1] - A[1])) % P)
    Sj = chord(acc, lad[j])
    if Sj is None: return None
    return ((T[0] - acc[0]) * (Sj[1] - acc[1]) - (T[1] - acc[1]) * (Sj[0] - acc[0])) % P, acc, Sj

def interp_roots(A, i, j, deg=8):
    """f(t1) is a rational function; sample it, interpolate the numerator, then root-find."""
    xs, ys = [], []
    t = 2
    while len(xs) < deg + 1 and t < 4000:
        r = resid(A, i, j, t)
        if r is not None: xs.append(t); ys.append(r[0])
        t += 1
    # Lagrange -> coefficients (dense, tiny)
    n = len(xs); C = [0] * n
    for k in range(n):
        num = [1]; den = 1
        for m in range(n):
            if m == k: continue
            num = [0] + num[:] if False else [(num[q - 1] if q else 0) - xs[m] * (num[q] if q < len(num) else 0) for q in range(len(num) + 1)]
            den = den * (xs[k] - xs[m]) % P
        inv = pow(den, P - 2, P)
        for q in range(len(num)): C[q] = (C[q] + ys[k] * num[q] % P * inv) % P if q < n else C[q]
    while len(C) > 1 and C[-1] % P == 0: C.pop()
    return C

def polyroots(C):
    """roots of a low-degree poly over GF(P) by gcd(t^P - t, f) then trial split."""
    C = [c % P for c in C]
    while len(C) > 1 and C[-1] == 0: C.pop()
    d = len(C) - 1
    if d <= 0: return []
    def pmod(a, b):
        a = a[:]
        while len(a) >= len(b) and any(a):
            if a[-1] == 0: a.pop(); continue
            f = a[-1] * pow(b[-1], P - 2, P) % P; s = len(a) - len(b)
            for k in range(len(b)): a[s + k] = (a[s + k] - f * b[k]) % P
            while len(a) > 1 and a[-1] == 0: a.pop()
        return a
    def pmul(a, b):
        r = [0] * (len(a) + len(b) - 1)
        for x, u in enumerate(a):
            if u:
                for y, w in enumerate(b): r[x + y] = (r[x + y] + u * w) % P
        return r
    # t^P mod f  by square-and-multiply
    f = C; r = [0, 1]; res = [1]; e = P
    while e:
        if e & 1: res = pmod(pmul(res, r), f)
        r = pmod(pmul(r, r), f); e >>= 1
    g = [(res[k] if k < len(res) else 0) - (1 if k == 1 else 0) for k in range(max(len(res), 2))]
    g = [x % P for x in g]
    def pgcd(a, b):
        while any(b[:-1]) or (b and b[-1]):
            a, b = b, pmod(a, b)
            if len(b) == 1: break
        return a
    h = pgcd(f[:], g)
    return h

for name, A in (('seed=L0', lad[0]),):
    for (i, j, flo) in ((8, 132, 39029), (73, 132, 39027), (8, 73, 39027), (132, 218, 39027)):
        C = interp_roots(A, i, j)
        h = polyroots(C)
        print('%-8s leaves(%3d,%3d) floor %d : elimination degree %d ; gcd(t^P-t,f) degree %d -> %s'
              % (name, i, j, flo, len(C) - 1, len(h) - 1,
                 'ROOTS EXIST' if len(h) - 1 >= 1 else 'no root'), flush=True)
