"""TASK 1c: make the frame-B budget EXHAUSTIVE, not budgeted.

Two facts, both computed exactly here:
 (1) The 168 satisfied rows are HOMOGENEOUS (rhs 0), so the admissible knob deltas keeping a
     kept-set KEEP satisfied are exactly ker_Z(A_KEEP) = Z^34 cap ker_Q(A_KEEP).  Therefore if
     rank_Q(A_KEEP) == rank_Q(A_SAT) the integer feasible set is IDENTICAL: breaking those rows
     buys nothing at any budget.
 (2) rank(A_SAT) = 26 = rank(A_R162) + |ESS| = 20 + 6, and ESS = {2554,6816,8124,9123,9421,S}
     are the only rows whose single deletion drops the rank.

PACKING LEMMA (proved by explicit witness below): if there are t pairwise-disjoint subsets
D_1..D_t of the 162 redundant rows, each of rank 20, then for any B = Bess u Bred with
|Bred| < t we have SAT\B contains (ESS\Bess) u D_i for some i, whose rank is 20 + |ESS\Bess|
= 26 - |Bess| = rank(A_SAT\Bess).  So ker_Q(A_SAT\B) = ker_Q(A_SAT\Bess), and the integer
feasible set depends ONLY on Bess.

Consequence: for every budget b < t, minbreak(P) = min{|Bess| : Bess subset of ESS,
(SAT\Bess) + P integer-feasible}.  That is 2^6 = 64 lattices, not C(168,b) subsets.
"""
import sys, os, itertools, json, time
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup as S

KN = S.KNOB
ESS = [2554, 6816, 8124, 9123, 9421, 'S']
R162 = [e for e in S.SAT if e not in ESS]
assert len(R162) == 162


def vec(e):
    return [Fraction(S.rows[e].get(u, 0)) for u in KN]


def rank_of(rows):
    basis = []
    for v in rows:
        w = v[:]
        for bp, bv in basis:
            if w[bp] != 0:
                f = w[bp]; w = [a - f * b for a, b in zip(w, bv)]
        p = next((i for i, x in enumerate(w) if x != 0), None)
        if p is not None:
            d = w[p]; basis.append((p, [x / d for x in w]))
    return len(basis)


def extract_rank20(pool, want):
    """greedily pull a minimal subset of `pool` reaching rank `want`; returns (subset, rest)."""
    basis, taken, rest = [], [], []
    for e in pool:
        if len(taken) >= want:
            rest.append(e); continue
        w = vec(e)
        for bp, bv in basis:
            if w[bp] != 0:
                f = w[bp]; w = [a - f * b for a, b in zip(w, bv)]
        p = next((i for i, x in enumerate(w) if x != 0), None)
        if p is None:
            rest.append(e); continue
        d = w[p]; basis.append((p, [x / d for x in w])); taken.append(e)
    return taken, rest, len(basis)


t0 = time.time()
rS = rank_of([vec(e) for e in S.SAT]); rR = rank_of([vec(e) for e in R162])
rE = rank_of([vec(e) for e in ESS])
print('rank(A_SAT)=%d  rank(A_R162)=%d  rank(A_ESS)=%d  -> 20+6=%d' % (rS, rR, rE, rR + rE))
assert rS == 26 and rR == 20 and rE == 6

pool = list(R162); packs = []
while True:
    taken, rest, r = extract_rank20(pool, 20)
    if r < 20: break
    packs.append(taken); pool = rest
print('disjoint rank-20 subsets of the 162 redundant rows: t = %d  (sizes %s, %d rows unused)'
      % (len(packs), [len(p) for p in packs], len(pool)))
for i, p in enumerate(packs):
    assert rank_of([vec(e) for e in p]) == 20
    assert rank_of([vec(e) for e in p + ESS]) == 26
print('all %d packs verified: rank 20 alone, rank 26 with ESS' % len(packs))
T = len(packs)
print('=> for every budget b <= %d, breaking redundant rows is PROVABLY WORTHLESS.' % (T - 1))

# now the complete answer over the 64 essential-row lattices
print('\nenumerating all 2^6 = 64 essential break-sets x all 127 bought-sets ...')
res = {}
best = (0, None)
for k in range(1, 8):
    for P in itertools.combinations(S.FAIL, k):
        mb, mbB = None, None
        for nb in range(0, 7):
            if mb is not None: break
            for Bess in itertools.combinations(ESS, nb):
                keep = [e for e in S.SAT if e not in Bess]
                if S.solve(keep + list(P)) is not None:
                    mb, mbB = nb, list(Bess); break
        gain = (k - mb) if mb is not None else None
        res[str(list(P))] = {'minbreak': mb, 'B': [str(x) for x in mbB] if mbB else None, 'gain': gain}
        if gain is not None and gain > best[0]: best = (gain, list(P), mbB)
        print('  buy %-34s minbreak=%s  B=%s  GAIN=%s'
              % (list(P), mb, [str(x) for x in mbB] if mbB else None, gain), flush=True)
json.dump({'t_packs': T, 'ESS': [str(x) for x in ESS], 'res': res}, open('w_exhaust.json', 'w'), indent=1)
print('\nBEST GAIN over all 127 bought-sets and all 64 essential break-sets: %s' % (best,))
print('elapsed %.0fs' % (time.time() - t0))
