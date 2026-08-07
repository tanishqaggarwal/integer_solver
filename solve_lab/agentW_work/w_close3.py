"""W round 2, item 1: THE CLOSING TEST on K+ = 40, from the EXACT cocircuit enumeration.

Lemma (round 1, restated for 40 knobs / 198 homogeneous SAT rows -- homogeneity re-verified in
w_hom2.py).  Any knob delta x breaks a set T = N(x) which is a UNION OF MINIMAL COCIRCUITS, so
testing unions of minimal cocircuits of total size <= 6 covers every j <= 7.

Two reductions used here, both exact:
 (R1) ker_Q(A_SAT\B) is monotone in B, and integer feasibility of (SAT\B)+P is monotone in the
      kernel.  So for a bought-set of size k it suffices to test the INCLUSION-MAXIMAL break
      sets of size <= k-1.
 (R2) every minimal support from the mod-2^61-1 search is re-verified to be genuinely
      rank-dropping over Q; mod-p artefacts are dropped and counted.
"""
import sys, os, itertools, json, time
from fractions import Fraction
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup2 as S

raw = []; meta = []
FILES = os.environ.get('WFILES', 'w_cocirc3_raw_s1_3.json,w_cocirc3_raw_s4_6.json').split(',')
for f in FILES:
    d = json.load(open(f)); raw += d['supports']; meta.append((f, len(d['supports']), d['subspace_hits']))
    print('%s: %d supports (positive-dimensional nodes %d)' % (f, len(d['supports']), d['subspace_hits']))
cc = {frozenset(int(x) if x != 'S' else 'S' for x in T) for T in raw}
# subset-minimality in O(n * 2^6): c is minimal iff no PROPER subset of it was also recorded
def _minimal(sets):
    out = []
    for c in sets:
        L = sorted(c, key=str)
        if any(frozenset(z) in sets for r in range(len(L)) for z in itertools.combinations(L, r)):
            continue
        out.append(c)
    return out
minimal = sorted(_minimal(cc), key=lambda c: (len(c), sorted(map(str, c))))
print('candidate supports %d ; subset-minimal %d ; by size %s'
      % (len(cc), len(minimal), dict(sorted(Counter(len(c) for c in minimal).items()))), flush=True)

KN = S.KNOB
V = lambda e: [Fraction(S.rows[e].get(u, 0)) for u in KN]
VEC = {e: V(e) for e in S.SAT}
def reduce_into(basis, v):
    w = v[:]
    for bp, bv in basis:
        if w[bp] != 0:
            f = w[bp]; w = [a-f*b for a, b in zip(w, bv)]
    p = next((i for i, x in enumerate(w) if x != 0), None)
    if p is None: return False
    dd = w[p]; basis.append((p, [x/dd for x in w])); return True
def rank_of(rows):
    basis = []
    for v in rows: reduce_into(basis, v)
    return len(basis)
rS = rank_of([VEC[e] for e in S.SAT])
# a spanning set of 32 independent rows, used to reject non-rank-dropping supports fast
SPAN, _b = [], []
for e in S.SAT:
    if reduce_into(_b, VEC[e]): SPAN.append(e)
SPANSET = set(SPAN)
def drops_rank(c):
    if not (SPANSET & c): return False          # the whole basis survives -> rank unchanged
    b = []
    for e in SPAN:
        if e not in c: reduce_into(b, VEC[e])
    if len(b) == rS: return False
    for e in S.SAT:
        if e in c or e in SPANSET: continue
        if reduce_into(b, VEC[e]) and len(b) == rS: return False
    return len(b) < rS
print('rank_Q(A_SAT) = %d ; exact-filtering %d minimal supports ...' % (rS, len(minimal)), flush=True)
t0 = time.time(); genuine = []; bogus = 0
for c in minimal:
    if drops_rank(c): genuine.append(c)
    else: bogus += 1
print('  genuinely rank-dropping %d ; mod-p artefacts dropped %d  (%.0fs)'
      % (len(genuine), bogus, time.time()-t0), flush=True)
minimal = sorted(_minimal(set(genuine)), key=lambda c: (len(c), sorted(map(str, c))))
print('  MINIMAL COCIRCUITS: %d ; by size %s'
      % (len(minimal), dict(sorted(Counter(len(c) for c in minimal).items()))), flush=True)
for c in minimal:
    if len(c) <= 2: print('    size %d : %s' % (len(c), sorted(map(str, c))))

CAP = int(os.environ.get('WCAP', '2000000'))
Bs = {frozenset()}; frontier = {frozenset()}; capped = False
while frontier and not capped:
    nxt = set()
    for b in frontier:
        for c in minimal:
            u = b | c
            if len(u) <= 6 and u not in Bs:
                Bs.add(u); nxt.add(u)
                if len(Bs) > CAP: capped = True; break
        if capped: break
    frontier = nxt
print('\nunion closure (size <= 6): %d break-sets%s ; by size %s'
      % (len(Bs), '  [CAPPED]' if capped else '', dict(sorted(Counter(len(b) for b in Bs).items()))), flush=True)

bysize = {}
for b in Bs: bysize.setdefault(len(b), []).append(b)
def maximal_upto(n):
    pool = [b for sz in range(n+1) for b in bysize.get(sz, [])]
    bigs = [b for sz in range(n, -1, -1) for b in bysize.get(sz, [])]
    out = []
    for b in pool:
        if not any(b < d for d in bigs): out.append(b)
    return out
t0 = time.time(); res = {}; best = (0, None, None); n = 0
for k in range(1, 8):
    cand = maximal_upto(k-1)
    mb, mbB = None, None
    for Pb in itertools.combinations(S.FAIL, k):
        mb, mbB = None, None
        for b in sorted(cand, key=len):
            n += 1
            if S.solve([e for e in S.SAT if e not in b] + list(Pb)) is not None:
                mb, mbB = len(b), sorted(b, key=str); break
        gain = (k-mb) if mb is not None else None
        res[str(list(Pb))] = {'minbreak_lt_k': mb, 'B': [str(x) for x in mbB] if mbB else None, 'gain': gain}
        if gain and gain > best[0]:
            best = (gain, list(Pb), mbB)
            S.price(S.solve([e for e in S.SAT if e not in mbB] + list(Pb)),
                    'buy %s break %s' % (list(Pb), mbB), tagfile='close3')
        print('  buy %-42s (|maximal break-sets tested| = %d) -> %s'
              % (list(Pb), len(cand), 'GAIN %d B=%s' % (gain, [str(x) for x in mbB]) if gain else 'NO'), flush=True)
print('\n%d integer solves in %.0fs' % (n, time.time()-t0))
print('BEST GAIN over every bought-set and every rank-dropping break-set of size <= 6: %s' % (best,))
json.dump({'files': FILES, 'nminimal': len(minimal), 'nBs': len(Bs), 'capped': capped, 'bogus': bogus,
           'best_gain': best[0], 'res': res}, open(os.environ.get('WOUT','w_close3.json'), 'w'), indent=1)
