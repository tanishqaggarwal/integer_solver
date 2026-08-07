#!/usr/bin/env python3
"""Awkward-case planted answers for the m<=8 run (|alpha|=4 / |beta|=4, the ONLY split
that reaches m=8 with an a<=4 table and a b=4 scan).

The four cases are chosen to be exactly the ones that would break if the bookkeeping were
subtly wrong -- including the one that already caught my own predictor:

  allneg   every sign negative, so the table part's lowest-exponent term is negative and the
           global-negation WLOG has to be absorbed by the +- branch of the x-only match.
  edges    exponents 0 and 255 both present: the extremes of the ladder.
  adjacent adjacent exponent pairs (e, e+1) inside BOTH halves: stresses the strictly-
           increasing `nxt(s) = ((s>>1)+1)<<1` step, which is what forbids repeated exponents.
  tight    all eight exponents inside a narrow high window: stresses the recursion's
           `nxt(s) > 512-2*(SZ-1-depth)` pruning, the one place a too-eager cut loses answers.

For each plant we predict the exact engine line and the s0 it must appear at, so validation
costs one s0 slice instead of a full b=4 sweep.
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'aa_offsets.json')))
p = int(D['p']); N = int(D['N']); A_ = int(D['a'])
G = (int(D['G'][0]), int(D['G'][1]))
lad = [(int(a), int(b)) for a, b in
       json.load(open(os.path.join(HERE, '..', 'agentX_work', 'xdata.json')))['ladder']]
def inv(z): return pow(z, p - 2, p)
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if (x1 - x2) % p == 0:
        if (y1 + y2) % p == 0: return None
        l = (3 * x1 * x1 + A_) % p * inv(2 * y1 % p) % p
    else: l = (y2 - y1) % p * inv((x2 - x1) % p) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)
def mul(k, P):
    k %= N; R = None; Q = P
    while k > 0:
        if k & 1: R = add(R, Q)
        Q = add(Q, Q); k >>= 1
    return R
def sp(s):
    P = lad[s >> 1]
    return P if not (s & 1) else (P[0], (-P[1]) % p)
def idx(e, sg): return 2 * e + (0 if sg > 0 else 1)
M64 = (1 << 64) - 1
LAD_TXT = ''.join('%d %d\n' % (x, y) for x, y in lad)

CASES = {
 'allneg':  [(11,-1),(37,-1),(66,-1),(94,-1),(131,-1),(168,-1),(205,-1),(244,-1)],
 'edges':   [(0,1),(3,-1),(29,1),(77,-1),(150,1),(201,-1),(254,1),(255,-1)],
 'adjacent':[(40,1),(41,-1),(88,-1),(89,1),(160,1),(161,1),(230,-1),(231,-1)],
 'tight':   [(238,1),(239,-1),(241,1),(244,-1),(247,1),(250,1),(252,-1),(255,1)],
}
os.makedirs(os.path.join(HERE, 'data'), exist_ok=True)
out = []
for name, terms in CASES.items():
    terms = sorted(terms)
    k = sum(s << e for e, s in terms) % N
    Tp = mul(k, G)
    base = Tp                                   # offset c = 0 for these
    assert base == mul(k, G)
    alpha, beta = terms[:4], terms[4:]
    V = (1 if alpha[0][1] > 0 else -1) * sum(s << e for e, s in alpha)
    bidx = sorted(idx(e, -s) for e, s in beta)  # scan side = negated complementary terms
    Q = base
    for s in bidx: Q = add(Q, sp(s))
    code = 0
    for i, s in enumerate(bidx): code |= s << (16 * i)
    line = (4, code, bidx[-1], Q[0] & M64)
    s0 = bidx[0]
    # decode check, in Python, before the engine ever runs
    Sb = sum((-1 if (s & 1) else 1) << (s >> 1) for s in bidx)
    dec = [br for br in (1, -1) if (br * V - Sb) % N == k % N]
    assert dec, 'self-consistency failed for %s' % name
    with open(os.path.join(HERE, 'data', 'q_%s.txt' % name), 'w') as f:
        f.write('%d %d\n' % base); f.write(LAD_TXT)
    out.append({'name': name, 'terms': terms, 'k': str(k), 'unsigned_wt': bin(k).count('1'),
                'expect': {'sz': 4, 'code': code, 'slast': bidx[-1], 'key': Q[0] & M64},
                's0': s0, 'branch': dec[0]})
    print('%-9s m=8 lowest-sign=%+d  unsigned wt=%3d  s0=%3d  branch=%+d  key=%d'
          % (name, terms[0][1], bin(k).count('1'), s0, dec[0], Q[0] & M64))
json.dump(out, open(os.path.join(HERE, 'plants8.json'), 'w'), indent=1)
print('\nwrote 4 awkward m=8 plants -> data/q_*.txt, plants8.json')
