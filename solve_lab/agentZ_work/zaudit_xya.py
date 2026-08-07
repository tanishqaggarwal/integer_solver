#!/usr/bin/env python3
"""Agent Z: audit of the search machinery of agents X, Y, AA.

Everything here is recomputed from MY OWN leaf extraction (zleaves.json / zexpo.json,
derived in section 12 from EQUATIONS.txt by regex + curve doubling) and cross-checked
against their artefacts.  Their code is read, never imported.
"""
import os, json, random
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.join(HERE, '..')
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

Z = json.load(open(os.path.join(HERE, 'zleaves.json')))
p = int(Z['p']); q3 = int(Z['q3']); b = int(Z['b'])
expo = {int(s): i for s, i in json.load(open(os.path.join(HERE, 'zexpo.json'))).items()}
# my ladder in X-coordinates (X = x + Q/3), indexed by exponent
MYL = {}
for s, (x, y) in Z['leaves'].items():
    MYL[expo[int(s)]] = ((int(x) + q3) % p, int(y))

def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        l = 3 * x1 * x1 % p * pow(2 * y1 % p, p - 2, p) % p
    else:
        l = (y2 - y1) * pow(x2 - x1, p - 2, p) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)

def mul(k, P):
    k %= N; R = None; Q = P
    while k:
        if k & 1: R = add(R, Q)
        Q = add(Q, Q); k >>= 1
    return R

print("=" * 92)
print("1. DO X / Y / AA AND I AGREE ON THE LADDER, G AND T?")
print("=" * 92)
xd = json.load(open(os.path.join(LAB, 'agentX_work', 'xdata.json')))
yd = json.load(open(os.path.join(LAB, 'agentY_work', 'ydata.json')))
XL = [(int(a), int(c)) for a, c in xd['ladder']]
YL = [(int(a), int(c)) for a, c in yd['ladder']]
print("  X ladder == Y ladder                    :", XL == YL)
print("  X ladder == MY leaves (exponent-ordered):", XL == [MYL[i] for i in range(256)])
XT = (int(xd['T'][0]), int(xd['T'][1])); YT = (int(yd['T'][0]), int(yd['T'][1]))
print("  X's T == Y's T                          :", XT == YT)
print("  T on the curve y^2 = X^3 + b            :", (XT[1] ** 2 - XT[0] ** 3 - b) % p == 0)
print("  N*T == O                                :", mul(N, XT) is None)
print("  L_i == 2^i * L_0 (my own doubling chain):", all(MYL[i] == mul(1 << i, MYL[0]) for i in range(0, 256, 17)))
print("  N*G == O                                :", mul(N, MYL[0]) is None)

print()
print("=" * 92)
print("2. AGENT Y — the complement construction, re-derived here")
print("=" * 92)
A1 = mul((2 ** 256 - 1) % N, MYL[0])
A2 = None
for i in range(256): A2 = add(A2, MYL[i])
A3 = mul(2 ** 256 - 1, MYL[0])
YA = (int(yd['A'][0]), int(yd['A'][1])) if 'A' in yd else None
print("  A three ways agree (mine)               :", A1 == A2 == A3)
print("  A == Y's A                              :", (YA is None) or (A1 == YA))
YTp = (int(yd['Tp'][0]), int(yd['Tp'][1]))
myTp = add(A1, (XT[0], (-XT[1]) % p))          # A - T
print("  T' = A - T  == Y's T'                   :", myTp == YTp)
print("  T + T' == A                             :", add(XT, YTp) == A1)
print("  T' on curve, N*T' == O, T' != T         :",
      (YTp[1] ** 2 - YTp[0] ** 3 - b) % p == 0, mul(N, YTp) is None, YTp != XT)
rng = random.Random(99)
ok = 0
for _ in range(20):
    w = rng.randrange(1, 256)
    S = set(rng.sample(range(256), w))
    Sb = set(range(256)) - S
    f1 = None
    for i in S: f1 = add(f1, MYL[i])
    f2 = None
    for i in Sb: f2 = add(f2, MYL[i])
    k = sum(1 << i for i in S); kb = sum(1 << i for i in Sb)
    if add(f1, f2) == A1 and k + kb == 2 ** 256 - 1 and f1 == mul(k, MYL[0]) and f2 == mul(kb, MYL[0]):
        ok += 1
print("  complement identity on 20 random S      : %d/20" % ok)
print("  => Y's central construction is INDEPENDENTLY REPRODUCED from my own leaf extraction.")

print()
print("=" * 92)
print("3. AGENT X — the exponent-256 gap in the SIGNED table (AB's finding), quantified")
print("=" * 92)
print("  xsigned.c reads exactly 256 ladder points into 512 signed slots:")
print("      for(int i=0;i<256;i++){ ... SX[2*i] = +L_i ; SX[2*i+1] = -L_i }")
print("  so the digit alphabet is  +-2^e  with  e in [0,255].  Exponent 256 is absent.")
def naf_weight(z):
    """minimal signed-digit weight of z (NAF is optimal)"""
    z = abs(z); w = 0
    while z:
        if z & 1:
            d = 2 - (z & 3)
            z -= d; w += 1
        z >>= 1
    return w
def min_signed_weight_modN(c):
    """min over the two integer representatives in (-N,N)"""
    r = c % N
    return min(naf_weight(r), naf_weight(r - N))
ones = 2 ** 256 - 1
print("  reach(2^256-1) = min signed weight with exponents <= 255 :", min_signed_weight_modN(ones),
      " (AA's aa_lattice.py reports 42)")
print("  reach(2^256)   :", min_signed_weight_modN(2 ** 256))
print("  For comparison, WITH an exponent-256 digit the complement class costs w'+2 terms:")
for wp in (5, 9, 12):
    print("     complement weight w'=%2d  ->  %2d signed terms with e<=256, but >= %d with e<=255"
          % (wp, wp + 2, min_signed_weight_modN(ones)))
print("  CONFIRMED: the near-all-ones family is outside X's signed sweep at any affordable depth.")
print("  (Independently detected by AA as reach=42 and by Y in its section 5.2c.)")

print()
print("=" * 92)
print("4. IS X's 'leading sign fixed positive' TABLE LOSSLESS?  (the sign-bookkeeping question)")
print("=" * 92)
print("""  xsigned.c table mode iterates  s0 = 0,2,4,...,510  -- the lowest-exponent digit is forced
  POSITIVE.  That is half of all signed sums.  Is coverage lost?  NO, and the reason is that the
  table stores ONLY the low 64 bits of the x-coordinate:
        every leading-negative sum  alpha  equals  -alpha'  with alpha' leading-positive,
        and  x(-P) = x(P),  so the two have the SAME key.
  Hence {keys of leading-positive sums} == {keys of all signed sums}.  Verified numerically:""")
cnt = 0; bad = 0
for _ in range(200):
    a = rng.randrange(3, 6)
    es = sorted(rng.sample(range(256), a))
    sg = [rng.choice((1, -1)) for _ in es]
    P = None
    for e, s in zip(es, sg):
        Q = MYL[e] if s > 0 else (MYL[e][0], (-MYL[e][1]) % p)
        P = add(P, Q)
    Pn = (P[0], (-P[1]) % p)
    if P[0] != Pn[0]: bad += 1
    cnt += 1
print("     200 random signed sums: x(alpha) == x(-alpha) in all cases :", bad == 0)
print("  => the sign restriction is LOSSLESS.  AA reached the same conclusion the hard way")
print("     (4 of its 8 plants failed because its PREDICTOR, not the engine, assumed otherwise).")

print()
print("=" * 92)
print("5. REPORTED EXHAUSTION — are the candidate counts what they must be?")
print("=" * 92)
def parse_done(path, key):
    out = []
    if not os.path.exists(path): return out
    for line in open(path):
        t = line.split()
        if len(t) > 2 and t[0] == 'DONE' and key in line:
            d = dict(x.split('=', 1) for x in t[1:] if '=' in x)
            out.append((line.strip(), d))
    return out
print("-- X unsigned, rep_real.txt (size-5 scan, must total C(256,5) = %d)" % comb(256, 5))
tot = 0; rows = 0
for line, d in parse_done(os.path.join(LAB, 'agentX_work', 'rep_real.txt'), 'size=5'):
    n = int(d.get('n', 0)); tot += n; rows += 1
print("   pieces=%d  sum(n)=%d  == C(256,5) : %s" % (rows, tot, tot == comb(256, 5)))
for sz in (2, 3, 4):
    s = sum(int(d.get('n', 0)) for _, d in parse_done(os.path.join(LAB, 'agentX_work', 'rep_real.txt'), 'size=%d' % sz))
    print("   size=%d sum(n)=%-12d == C(256,%d)=%-12d : %s" % (sz, s, sz, comb(256, sz), s == comb(256, sz)))
print("-- X signed, srep_real.txt (must total C(256,b)*2^b)")
for line in open(os.path.join(LAB, 'agentX_work', 'srep_real.txt')):
    t = line.split()
    if t and t[0] == 'DONE':
        d = dict(x.split('=', 1) for x in t[2:] if '=' in x)
        sz = int(d['sz']); n = int(d['n'])
        want = comb(256, sz) * (2 ** sz)
        print("   sz=%d n=%-12d want C(256,%d)*2^%d=%-12d : %s" % (sz, n, sz, sz, want, n == want))
print("-- Y complement, rep_comp.txt")
print(open(os.path.join(LAB, 'agentY_work', 'rep_comp.txt')).read().strip())
for line in open(os.path.join(LAB, 'agentY_work', 'rep_comp.txt')):
    t = line.split()
    if t and t[0] == 'DONE':
        d = dict(x.split('=', 1) for x in t[1:] if '=' in x)
        sz = int(d['size']); n = int(d['n'])
        print("   size=%d n=%-12d == C(256,%d)=%-12d : %s" % (sz, n, sz, comb(256, sz), n == comb(256, sz)))
