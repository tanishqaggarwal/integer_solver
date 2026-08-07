"""Extend the §79 superset result from the complete 2^12 to the complete 2^16.

§79 proved, exhaustively over 2^12: every subset attaining 39,026 CONTAINS the witness
{642, 28730, 29854, 31864}, but containing it does not guarantee 39,026 (7 of 28 supersets
at |W|=6 fall below).  The 2^16 space is four handles wider, so the claim has to be re-tested
there rather than assumed -- the extra handles could in principle carry an optimum of their
own.

This re-prices only the subsets that need it: every superset of the witness inside H16, and
(as the converse test) every subset that scored 39,026 according to the completed checkpoint.
"""
import sys, os, pickle, itertools, collections, json

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
os.chdir(MDIR); sys.path.insert(0, MDIR)
import shim                                                    # noqa: F401
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
import ieng                                                    # noqa: E402

W4 = (642, 28730, 29854, 31864)
PF = json.load(open('pfamily.json'))
H16 = sorted({v['h'] for v in PF['incident_12'].values()})
ck = pickle.load(open('enumsub16.pkl', 'rb'))
assert ck['complete'], 'checkpoint is not complete'
print(f'H16 = {H16}  ({len(H16)} handles), checkpoint complete at {ck["i"]:,}')

ORDER = [W for k in range(len(H16) + 1) for W in itertools.combinations(H16, k)]
sup = [W for W in ORDER if set(W4) <= set(W)]
print(f'supersets of the witness in H16: {len(sup):,}  (expect 2^{len(H16)-4} = {2**(len(H16)-4):,})')

winners = []
supwin = collections.Counter(); suptot = collections.Counter()
for W in sup:
    sc = ieng.tune(list(W))['score']
    suptot[len(W)] += 1
    if sc >= 39026:
        supwin[len(W)] += 1
        winners.append(W)
print(f'\nsupersets scoring >= 39026: {len(winners)} of {len(sup)}')
for k in sorted(suptot):
    print(f'  |W|={k:2d}: {supwin[k]:3d} winners of {suptot[k]:5d} supersets at that size')

# converse: are there winners that are NOT supersets?  the completed run counted 114 at 39,026,
# so if the supersets account for all 114 the containment claim holds over 2^16.
n114 = sum(c.get(39026, 0) for c in ck['bysize'].values())
print(f'\ncheckpoint total at 39,026 : {n114}')
print(f'supersets attaining it     : {len(winners)}')
print(f'therefore NON-supersets    : {n114 - len(winners)}')
print('CLAIM VERIFIED over 2^16' if n114 == len(winners) else 'CLAIM FAILS over 2^16')

json.dump({'space': '2^16', 'handles': H16, 'witness': list(W4),
           'n_supersets': len(sup), 'n_winners': len(winners),
           'checkpoint_total_39026': n114,
           'non_supersets_attaining': n114 - len(winners),
           'by_size': {str(k): [supwin[k], suptot[k]] for k in sorted(suptot)},
           'winners': [list(w) for w in winners]},
          open('verifysup16.json', 'w'), indent=1)
