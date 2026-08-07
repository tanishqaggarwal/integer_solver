#!/usr/bin/env python3
"""agent AE -- INDEPENDENT verifier for any candidate scalar.

Deliberately shares no code with ae_lib.py:
  * Jacobian coordinates, not affine        (different formulas)
  * fixed-window (w=4) scalar multiplication, not binary double-and-add
  * curve parameters re-read from agentX_work/xdata.json, not from ae_data.json
  * modular inverse by extended Euclid, not by Fermat's little theorem
so an error in one implementation cannot hide in the other.

  python3 ae_verify.py <decimal scalar> [...]        verify k*G == T, print popcount
  python3 ae_verify.py --selftest                    prove the verifier on known values
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
_d = json.load(open(os.path.join(HERE, '..', 'agentX_work', 'xdata.json')))
P = int(_d['p']); B = int(_d['b']); N = int(_d['N'])
Gx, Gy = int(_d['G'][0]), int(_d['G'][1])
Tx, Ty = int(_d['T'][0]), int(_d['T'][1])

def egcd_inv(a, m):
    """extended Euclid -- not pow(a, m-2, m)"""
    a %= m
    old_r, r = a, m
    old_s, s = 1, 0
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
    if old_r != 1: raise ZeroDivisionError('not invertible')
    return old_s % m

# ---- Jacobian arithmetic on y^2 = x^3 + B  (a = 0) -------------------------
INF = (0, 1, 0)

def j_dbl(Pt):
    X, Y, Z = Pt
    if Z == 0 or Y == 0: return INF
    A = (4 * X * Y * Y) % P
    Bq = (8 * pow(Y, 4, P)) % P
    C = (3 * X * X) % P                       # a = 0
    X3 = (C * C - 2 * A) % P
    Y3 = (C * (A - X3) - Bq) % P
    Z3 = (2 * Y * Z) % P
    return (X3, Y3, Z3)

def j_add(Pt, Qt):
    X1, Y1, Z1 = Pt; X2, Y2, Z2 = Qt
    if Z1 == 0: return Qt
    if Z2 == 0: return Pt
    Z1Z1 = Z1 * Z1 % P; Z2Z2 = Z2 * Z2 % P
    U1 = X1 * Z2Z2 % P; U2 = X2 * Z1Z1 % P
    S1 = Y1 * Z2 % P * Z2Z2 % P; S2 = Y2 * Z1 % P * Z1Z1 % P
    if U1 == U2:
        if S1 != S2: return INF
        return j_dbl(Pt)
    H = (U2 - U1) % P
    I = (2 * H) ** 2 % P
    J = H * I % P
    r = 2 * (S2 - S1) % P
    V = U1 * I % P
    X3 = (r * r - J - 2 * V) % P
    Y3 = (r * (V - X3) - 2 * S1 * J) % P
    Z3 = (((Z1 + Z2) ** 2 - Z1Z1 - Z2Z2) * H) % P
    return (X3, Y3, Z3)

def j_affine(Pt):
    X, Y, Z = Pt
    if Z == 0: return None
    zi = egcd_inv(Z, P); zi2 = zi * zi % P
    return (X * zi2 % P, Y * zi2 % P * zi % P)

def mul_window(k, Pt, w=4):
    """fixed-window scalar multiplication -- structurally different from binary double-and-add"""
    k %= N
    if k == 0: return INF
    tbl = [INF, Pt]
    for i in range(2, 1 << w):
        tbl.append(j_add(tbl[-1], Pt) if i & 1 else j_dbl(tbl[i >> 1]))
    bits = k.bit_length()
    nd = (bits + w - 1) // w
    acc = INF
    for i in range(nd - 1, -1, -1):
        if acc != INF:
            for _ in range(w): acc = j_dbl(acc)
        d = (k >> (i * w)) & ((1 << w) - 1)
        if d: acc = j_add(acc, tbl[d])
    return acc

def on_curve(x, y): return (y * y - x * x * x - B) % P == 0

def verify(k):
    Q = j_affine(mul_window(k, (Gx, Gy, 1)))
    ok = (Q == (Tx, Ty))
    return ok, Q

def selftest():
    global Tx, Ty
    fails = 0
    G = (Gx, Gy, 1)
    print('curve/param checks:')
    print('  G on curve             :', on_curve(Gx, Gy))
    print('  T on curve             :', on_curve(Tx, Ty))
    print('  N*G == infinity        :', mul_window(N, G)[2] % P == 0 or j_affine(mul_window(N, G)) is None)
    # 1..8 * G against a from-scratch repeated-addition chain in affine coords
    def aff_add(A, Bp):
        if A is None: return Bp
        if Bp is None: return A
        if A[0] == Bp[0]:
            if (A[1] + Bp[1]) % P == 0: return None
            l = 3 * A[0] * A[0] * egcd_inv(2 * A[1], P) % P
        else:
            l = (Bp[1] - A[1]) * egcd_inv(Bp[0] - A[0], P) % P
        x = (l * l - A[0] - Bp[0]) % P
        return (x, (l * (A[0] - x) - A[1]) % P)
    acc = None
    for i in range(1, 40):
        acc = aff_add(acc, (Gx, Gy))
        if j_affine(mul_window(i, G)) != acc: fails += 1
    print('  k*G matches a repeated-addition chain for k = 1..39 :', fails == 0)
    # ladder cross-check: 2^i * G for a spread of i
    lad = json.load(open(os.path.join(HERE, '..', 'agentX_work', 'xdata.json')))['ladder']
    bad = 0
    for i in (0, 1, 7, 31, 64, 127, 200, 255):
        if j_affine(mul_window(1 << i, G)) != (int(lad[i][0]), int(lad[i][1])): bad += 1
    print('  2^i*G matches the instance ladder for 8 spread i    :', bad == 0)
    fails += bad
    # end-to-end: plant a scalar, confirm verify() accepts only it
    kk = 0xDEADBEEFCAFE1234567890ABCDEF
    Qp = j_affine(mul_window(kk, G))
    sx, sy = Tx, Ty
    Tx, Ty = Qp
    good = verify(kk)[0] and not verify(kk + 1)[0]
    Tx, Ty = sx, sy
    print('  planted scalar accepted, k+1 rejected               :', good)
    fails += (0 if good else 1)
    print('SELFTEST', 'PASS' if fails == 0 else 'FAIL (%d)' % fails)
    return fails == 0

if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] == '--selftest':
        sys.exit(0 if selftest() else 1)
    for a in sys.argv[1:]:
        k = int(a)
        ok, Q = verify(k)
        print('k = %d' % k)
        print('  k*G == T                : %s' % ok)
        print('  popcount(k)             : %d' % bin(k).count('1'))
        print('  popcount(k mod N)       : %d' % bin(k % N).count('1'))
        print('  bit_length              : %d' % k.bit_length())
        if not ok: print('  computed k*G            : %r' % (Q,))
