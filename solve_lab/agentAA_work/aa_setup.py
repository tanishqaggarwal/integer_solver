#!/usr/bin/env python3
"""Agent AA -- offset-shifted signed-digit MITM.

  k*G = T   <=>   (k-c)*G = T - c*G

So a search for "k - c has signed-digit weight <= m" is the IDENTICAL search machinery
with the base point B replaced by B_c = T - c*G.  The MITM *table* (all signed a-term
ladder sums) is OFFSET-INDEPENDENT: build once, reuse for every offset.  Only the scan
side moves.  That asymmetry is what makes |C| offsets cost |C| x (cheap half) instead of
|C| x (whole search).

This script:
  1. re-verifies the curve / ladder / target from agentX_work/xdata.json (read-only)
     by independent recomputation (not by trusting the file),
  2. builds the offset list C with a computed redundancy measure for each,
  3. writes one engine data file per offset (base point = T - c*G),
  4. writes planted-answer data files for validation.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
X = os.path.join(HERE, '..', 'agentX_work')

d = json.load(open(os.path.join(X, 'xdata.json')))
p = int(d['p']); A_ = int(d['a']); B_ = int(d['b']); N = int(d['N'])
G = (int(d['G'][0]), int(d['G'][1]))
T = (int(d['T'][0]), int(d['T'][1]))
lad = [(int(x), int(y)) for x, y in d['ladder']]

def inv(z): return pow(z, p - 2, p)
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if (x1 - x2) % p == 0:
        if (y1 + y2) % p == 0: return None
        l = (3 * x1 * x1 + A_) % p * inv(2 * y1 % p) % p
    else:
        l = (y2 - y1) % p * inv((x2 - x1) % p) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)
def neg(P): return None if P is None else (P[0], (-P[1]) % p)
def mul(k, P):
    k %= N
    R = None; Q = P
    while k > 0:
        if k & 1: R = add(R, Q)
        Q = add(Q, Q); k >>= 1
    return R
def oncur(P): return P is None or (P[1] * P[1] - pow(P[0], 3, p) - A_ * P[0] - B_) % p == 0

# ---------------------------------------------------------------- 1. VERIFY
def verify():
    ok = True
    c1 = (p == 2**256 - 2**32 - 977)
    c2 = (N == 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141)
    print('p == 2^256-2^32-977          :', c1)
    print('N == secp256k1 group order   :', c2)
    ok &= c1 and c2
    # ladder independently recomputed by repeated doubling from G
    cur = G; bad = 0
    for i in range(256):
        if lad[i] != cur: bad += 1
        cur = add(cur, cur)
    print('L_i == 2^i*G (recomputed)    : %d/256 bad=%d' % (256 - bad, bad)); ok &= (bad == 0)
    onc = sum(1 for q in lad if oncur(q))
    print('ladder points on curve       : %d/256' % onc); ok &= (onc == 256)
    c3 = (mul(N, G) is None); print('N*G == O                     :', c3); ok &= c3
    c4 = oncur(T); print('T on curve                   :', c4); ok &= c4
    c5 = (mul(N, T) is None); print('N*T == O                     :', c5); ok &= c5
    # G is the exponent-0 ladder point
    c6 = (lad[0] == G); print('lad[0] == G                  :', c6); ok &= c6
    print('VERIFY:', 'PASS' if ok else 'FAIL')
    assert ok
verify()

# ---------------------------------------------------------------- 2. OFFSETS
def naf(n):
    """non-adjacent form: minimal-weight signed binary rep. returns list of (exp,sign)."""
    out = []; i = 0
    while n:
        if n & 1:
            z = 2 - (n & 3)
            out.append((i, 1 if z == 1 else -1)); n -= z
        n >>= 1; i += 1
    return out

def reach(c):
    """min number of signed terms with exponents in [0,255] needed to express c mod N
       (as +c or -c, over the two integer representatives in (-N,N)).  inf if it needs
       an exponent > 255.  This is exactly 'how many extra terms the plain c=0 search
       would have to spend to cover this offset for free'."""
    r = c % N
    best = float('inf')
    for v in (r, r - N):
        a = abs(v)
        if a == 0: return 0
        f = naf(a)
        if max(e for e, _ in f) <= 255:
            best = min(best, len(f))
    return best

M = (1 << 256) - 1
def rep(pat_bits, pat):        # repeating pattern
    v = 0; w = pat_bits
    for i in range(256 // w): v |= pat << (w * i)
    return v & M

OFF = []   # (tag, value, rationale, raw, tier)
def O(tag, val, why, tier=1): OFF.append((tag, val % N if val % N else 0, why, val, tier))

O('c0',        0,            'plain case: signed-digit weight of k itself (contains low Hamming weight AND low run-length)')
O('ones',      M,            'complement: k has few ZERO bits.  = all-ones constant. (agent Y owns this class)')
O('a55',       M // 3,       '0x5555.. = 01 repeating; the densest period-2 constant')
O('aAA',       2 * (M // 3), '0xAAAA.. = 10 repeating; the other period-2 phase')
O('a33',       M // 5,       '0x3333.. = 0011 repeating, period 4')
O('aCC',       4 * (M // 5), '0xCCCC.. = 1100 repeating, period 4 other phase')
O('a11',       M // 15,      '0x1111.. = 0001 repeating, period 4 sparse')
O('a0F',       M // 17,      '0x0F0F.. = 00001111 repeating, period 8')
O('a05',       M // 51,      '0x0505.. = period 8, two ones')
O('a03',       M // 85,      '0x0303.. = period 8, adjacent pair')
O('a01',       M // 255,     '0x0101.. = period 8, one bit per byte')
O('aFF00',     M // 257,     '0x00FF00FF.. = period 16, byte-alternating')
O('a0001',     M // 65535,   '0x0001.. = period 16, one bit per 16')
O('aFFFF',     M // 65537,   '0x0000FFFF.. = period 32')
O('a32_1',     M // (2**32 - 1),  'period-32, one bit per word (8 terms)')
O('a32_F',     M // (2**32 + 1),  'period-64 word-alternating')
O('a64_1',     M // (2**64 - 1),  'period-64, one bit per 64 (4 terms) -- near-redundant, kept as control')
O('halfN',     (N - 1) // 2, 'half the group order; "k near N/2" is the natural midpoint choice')
O('inv2',      pow(2, N - 2, N), '1/2 mod N -- T = G/2 style key')
O('inv3',      pow(3, N - 2, N), '1/3 mod N')
O('inv5',      pow(5, N - 2, N), '1/5 mod N')
O('inv7',      pow(7, N - 2, N), '1/7 mod N')
O('inv10',     pow(10, N - 2, N), '1/10 mod N')
O('lam',       0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72, 'GLV eigenvalue lambda: k = lambda + (small signed) is a+b*lambda with b=1 and a of low weight')
O('lam2',      (0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72**2) % N, 'lambda^2 mod N = the third cube root of 1')
O('pmodN',     p % N,        'the field prime read as a scalar (2^256-2^32-977 mod N)')
O('2p256',     (2**256) % N, '2^256 mod N: reaches k of the form 2^256 - (small), i.e. exponent 256 which the 256-point ladder cannot represent')
O('bcurve',    B_ % N,       'the instance curve constant b read as a scalar')
O('cshift',    int(json.load(open(os.path.join(X, '..', 'agentQ_work', 'curve.json')))['c_shift']) % N,
                             'the instance c_shift constant read as a scalar')
O('d1e76',     10**76,       'round DECIMAL constant 10^76 (~2^252): "nice" in base 10, structureless in base 2')
O('d1e38',     10**38,       'round decimal 10^38')
O('drep',      (10**77 - 1) // 9, 'decimal repunit 111...1 (77 digits)')

# ---- WAVE 2a: the exponent-256/257 terms the 256-point ladder structurally cannot hold.
# The basis loads exponents 0..255 only, so a key of the form 2^256 - (something small) has
# NO short signed representation in it.  Extending the basis is one way to fix that; an
# offset is a strictly better way, because offset c at depth m covers the e=256 branch with
# m further terms whereas a 257-point basis at depth m covers it with only m-1.
O('p2p256',    2**256,       'exponent 256: the term the 256-point ladder cannot represent at all', 1)
O('n2p256',    -(2**256),    '= 2N mod 2^256.  the NEGATIVE exponent-256 term (other sign branch)', 1)
O('p2p256p1',  2**256 + 1,   '2^256+1: round-above-the-top constant', 1)
O('n2p256m1',  -((2**256) - 1), 'the NEGATED all-ones constant (k = -(2^256-1) + small)', 1)
O('p2p257',    2**257,       'exponent 257', 2)
O('n2p257',    -(2**257),    'exponent 257, negative branch', 2)

# ---- WAVE 2b: negations of the structured constants.  c and -c are DIFFERENT offsets
# (offset c covers k = c + small; offset -c covers k = -c + small = N - c + small, i.e.
# "the group order minus a nice constant"), so each needs its own scan.
O('n_a55',     -(M // 3),    'N - 0x5555.. ', 2)
O('n_a33',     -(M // 5),    'N - 0x3333.. ', 2)
O('n_a11',     -(M // 15),   'N - 0x1111.. ', 2)
O('n_a0F',     -(M // 17),   'N - 0x0F0F.. ', 2)
O('n_a01',     -(M // 255),  'N - 0x0101.. ', 2)
O('n_aFF00',   -(M // 257),  'N - 0x00FF00FF.. ', 2)
O('n_a0001',   -(M // 65535),'N - 0x0001.. ', 2)
O('n_lam',     -0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72, 'N - lambda (= the other GLV root, -lambda-1 up to 1)', 2)
O('n_inv3',    -pow(3, N - 2, N), '-1/3 mod N', 2)
O('n_d1e76',   -(10**76),    'N - 10^76', 2)

O('N_noop',    N,            'CONTROL: c = N is identically 0 mod N -- must reproduce c0 exactly', 3)
O('p2_128',    2**128,       'CONTROL: single power of two -- provably redundant (reach=1)', 3)
O('p2_128p1',  2**128 + 1,   'CONTROL: 2^128+1 -- provably redundant (reach=2)', 3)
O('p2_255_1',  2**255 - 1,   'CONTROL: 2^255-1 -- provably redundant (reach=2)', 3)
O('p2_255',    2**255,       'CONTROL: 2^255 -- provably redundant (reach=1)', 3)

# de-duplicate on value mod N, keeping first tag
seen = {}
FINAL = []
for tag, v, why, raw, tier in OFF:
    if v in seen:
        print('DUP: %-10s == %s (same point mod N) -- dropped' % (tag, seen[v]))
        continue
    seen[v] = tag
    FINAL.append((tag, v, why, raw, tier))

print()
print('%-10s %-5s %-6s %-6s %s' % ('tag', 'tier', 'reach', 'nafwt', 'rationale'))
for tag, v, why, raw, tier in sorted(FINAL, key=lambda z: (z[4], -reach(z[1]))):
    r = reach(v)
    nw = len(naf(v % N)) if v % N else 0
    print('%-10s %-5d %-6s %-6d %s' % (tag, tier, ('inf' if r == float('inf') else r), nw, why[:88]))

# ---------------------------------------------------------------- 3. DATA FILES
LAD_TXT = ''.join('%d %d\n' % (x, y) for x, y in lad)
def write_data(path, base):
    with open(path, 'w') as f:
        f.write('%d %d\n' % base); f.write(LAD_TXT)

os.makedirs(os.path.join(HERE, 'data'), exist_ok=True)
manifest = []
for tag, v, why, raw, tier in FINAL:
    Bc = T if v % N == 0 else add(T, neg(mul(v, G)))
    assert Bc is not None and oncur(Bc)
    write_data(os.path.join(HERE, 'data', 'd_%s.txt' % tag), Bc)
    manifest.append({'tag': tag, 'c': str(v % N), 'tier': tier,
                     'reach': (None if reach(v) == float('inf') else reach(v)),
                     'nafwt': (len(naf(v % N)) if v % N else 0), 'why': why,
                     'base': [str(Bc[0]), str(Bc[1])]})

# sanity: c0 and N_noop must give the identical base point
b0 = [m for m in manifest if m['tag'] == 'c0'][0]['base']
bn = [m for m in manifest if m['tag'] == 'N_noop']
if bn: print('\nCONTROL c=N gives identical base to c=0 :', bn[0]['base'] == b0)

# ---------------------------------------------------------------- 4. PLANTS
import random
random.seed(20260807)
plants = []
PLANT_OFFSETS = ['c0', 'ones', 'a55', 'a01', 'inv3', 'lam', '2p256', 'd1e76',
                 'n2p256', 'n2p256m1', 'p2p256p1', 'p2p257', 'n2p257',
                 'n_a55', 'n_lam', 'n_a0001', 'halfN', 'a0F', 'aFFFF', 'drep',
                 'a64_1', 'p2_128', 'p2_255_1']
for tag in PLANT_OFFSETS:
    v = seen and [m for m in manifest if m['tag'] == tag][0]
    c = int(v['c'])
    exps = sorted(random.sample(range(256), 5))
    sgn = [random.choice((1, -1)) for _ in range(5)]
    delta = sum(s << e for e, s in zip(exps, sgn))
    kp = (c + delta) % N
    Tp = mul(kp, G)
    Bp = Tp if c % N == 0 else add(Tp, neg(mul(c, G)))
    # Bp must equal delta*G
    assert Bp == mul(delta % N, G), 'offset bookkeeping broken for %s' % tag
    write_data(os.path.join(HERE, 'data', 'p_%s.txt' % tag), Bp)
    plants.append({'tag': tag, 'c': str(c), 'exps': exps, 'sgn': sgn,
                   'delta': str(delta), 'k': str(kp),
                   'unsigned_wt': bin(kp).count('1'),
                   'base': [str(Bp[0]), str(Bp[1])]})
    print('plant %-7s m=5 delta terms=%s  k unsigned wt=%d' %
          (tag, list(zip(exps, sgn)), bin(kp).count('1')))

json.dump({'manifest': manifest, 'plants': plants,
           'T': [str(T[0]), str(T[1])], 'G': [str(G[0]), str(G[1])],
           'p': str(p), 'N': str(N), 'a': str(A_), 'b': str(B_)},
          open(os.path.join(HERE, 'aa_offsets.json'), 'w'), indent=1)
print('\nwrote %d real data files + %d plants -> data/, aa_offsets.json' % (len(manifest), len(plants)))
