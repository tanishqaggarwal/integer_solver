"""Verify a claim I made without proving it.

Round 13: "the other 29 subsets scoring 39,026 are all supersets of the witness".
The by-size counts were CONSISTENT with that but did not establish it -- a non-superset
could score 39,026 while some superset does not (and indeed at |W|=6 only 21 of the 28
supersets reached it, so 'all supersets win' is already false).

This records the actual winning subsets over the complete 2^12 space and checks each.
"""
import sys, os, json, time, itertools, collections, pickle
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import ieng

PF = json.load(open('pfamily.json'))
H12 = sorted({v['h'] for v in PF['incident_7'].values()})
WIT = {642, 28730, 29854, 31864}
print(f'H12 = {H12}', flush=True)
print(f'witness = {sorted(WIT)}', flush=True)

winners = []
n = 0
t0 = time.time()
for k in range(len(H12) + 1):
    for W in itertools.combinations(H12, k):
        n += 1
        try:
            r = ieng.tune(list(W)) if W else {'ok': True,
                                              'score': ieng.NEQ - len(ieng.FAILS_UNC)}
        except Exception:
            continue
        if r.get('ok') and r['score'] >= 39026:
            winners.append((W, r['score']))
    print(f'  |W|={k} done, {n} priced, winners so far {len(winners)} '
          f'({time.time()-t0:.0f}s)', flush=True)

print(f'\n{len(winners)} subsets score >= 39026 over the complete 2^12 space', flush=True)
sup = [W for W, s in winners if WIT <= set(W)]
non = [W for W, s in winners if not (WIT <= set(W))]
print(f'  supersets of the witness : {len(sup)}', flush=True)
print(f'  NOT supersets            : {len(non)}   {non}', flush=True)
print(f'\nCLAIM "all 39,026 optima are supersets of the witness": '
      f'{"VERIFIED" if not non else "FALSE"}', flush=True)

bysize = collections.Counter(len(W) for W, s in winners)
print('\nwinners by support size, against the number of supersets of that size:')
from math import comb
for k in sorted(bysize):
    nsup = comb(len(H12) - 4, k - 4) if k >= 4 else 0
    print(f'  |W|={k}: {bysize[k]} winners, {nsup} supersets of the witness exist at that size',
          flush=True)
json.dump({'winners': [list(W) for W, s in winners],
           'all_supersets': not non, 'non_supersets': [list(W) for W in non]},
          open('verifysup.json', 'w'), indent=1)
print('\nwrote verifysup.json', flush=True)
