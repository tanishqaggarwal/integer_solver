"""W round 2: re-run round-1's EXHAUSTIVE essential-row test on the EXTENDED knob set K+=40.
Essential rows are recomputed, not inherited."""
import sys, os, itertools, json, time
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup2 as S
KN = S.KNOB
def vec(e): return [Fraction(S.rows[e].get(u, 0)) for u in KN]
def rank_of(rows):
    basis = []
    for v in rows:
        w = v[:]
        for bp, bv in basis:
            if w[bp] != 0:
                f = w[bp]; w = [a - f*b for a, b in zip(w, bv)]
        p = next((i for i, x in enumerate(w) if x != 0), None)
        if p is not None:
            dd = w[p]; basis.append((p, [x/dd for x in w]))
    return len(basis)
t0 = time.time()
rS = rank_of([vec(e) for e in S.SAT])
ESS = [e for e in S.SAT if rank_of([vec(x) for x in S.SAT if x != e]) < rS]
RED = [e for e in S.SAT if e not in ESS]
print('knobs %d  SAT %d  rank(A_SAT)=%d' % (len(KN), len(S.SAT), rS))
print('ESSENTIAL rows (single deletion drops rank): %d  ->  %s' % (len(ESS), ESS))
print('redundant rows: %d   rank(A_RED)=%d' % (len(RED), rank_of([vec(e) for e in RED])))
# packing lemma: disjoint full-rank-on-RED subsets
def extract(pool, want):
    basis, taken, rest = [], [], []
    for e in pool:
        if len(taken) >= want: rest.append(e); continue
        w = vec(e)
        for bp, bv in basis:
            if w[bp] != 0:
                f = w[bp]; w = [a-f*b for a, b in zip(w, bv)]
        p = next((i for i, x in enumerate(w) if x != 0), None)
        if p is None: rest.append(e); continue
        dd = w[p]; basis.append((p, [x/dd for x in w])); taken.append(e)
    return taken, rest, len(basis)
rR = rank_of([vec(e) for e in RED]); pool = list(RED); packs = []
while True:
    tk, rest, r = extract(pool, rR)
    if r < rR: break
    packs.append(tk); pool = rest
print('disjoint rank-%d subsets of the %d redundant rows: t = %d  => redundant breaks are '
      'provably worthless for every budget b <= %d' % (rR, len(RED), len(packs), len(packs)-1))
print()
best = (0, None); res = {}
nE = len(ESS)
print('enumerating all 2^%d = %d essential break-sets x all 127 bought-sets ...' % (nE, 2**nE))
for k in range(1, 8):
    for P in itertools.combinations(S.FAIL, k):
        mb, mbB = None, None
        for nb in range(0, nE+1):
            if mb is not None: break
            for Bess in itertools.combinations(ESS, nb):
                keep = [e for e in S.SAT if e not in Bess]
                if S.solve(keep + list(P)) is not None:
                    mb, mbB = nb, list(Bess); break
        gain = (k-mb) if mb is not None else None
        res[str(list(P))] = {'minbreak': mb, 'B': [str(x) for x in mbB] if mbB else None, 'gain': gain}
        if gain is not None and gain > best[0]: best = (gain, list(P), mbB)
        print('  buy %-34s minbreak=%s B=%s GAIN=%s' % (list(P), mb, mbB, gain), flush=True)
print()
print('BEST GAIN over all 127 bought-sets and all %d essential break-sets: %s' % (2**nE, best,))
print('elapsed %.0fs' % (time.time()-t0))
json.dump({'KNOB': KN, 'ESS': [str(x) for x in ESS], 'packs': len(packs), 'res': res,
           'best_gain': best[0]}, open('w_exhaust2.json', 'w'), indent=1)
