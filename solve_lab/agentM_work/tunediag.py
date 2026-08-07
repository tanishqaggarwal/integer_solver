"""Why does the tuner fail calibration? Use the deliverable as a KNOWN-GOOD solution.

The deliverable's own handle values, expressed as a delta from the uncorrupted baseline,
are by construction a tuning that reaches 39,026.  If the greedy tuner cannot find it,
the question is whether that delta is even inside the space the tuner searches.
"""
import sys, os, collections, json
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as EB, engine3 as E3
import price as PR, fscore

vd = PR.load_deliverable()
P = PR.TunedPricer(vd)
h0 = [642, 28730, 29854, 31864]
freed, demote = PR.closure(h0)
eng = E3.Eng(demote)

# uncorrupted baseline inside this engine
seed_b = {f: vd[f] for f in eng.FREE if vd[f] != 0}
for u in freed:
    if P.v_unc[u]:
        seed_b[u] = P.v_unc[u]
    else:
        seed_b.pop(u, None)
v0 = eng.forward(seed_b)
bad0 = eng.badatoms(v0)
F0 = sorted(fscore.fails(bad0))
print(f'UNCORRUPTED BASELINE: score {fscore.score(bad0)}, {len(F0)} failures, {len(bad0)} bad atoms')
print(f'  bad atoms  : {sorted(bad0)}')
print(f'  failures   : {F0}')

# the deliverable, in the same engine
seed_d = {f: vd[f] for f in eng.FREE if vd[f] != 0}
vD = eng.forward(seed_d)
badD = eng.badatoms(vD)
FD = sorted(fscore.fails(badD))
print(f'\nDELIVERABLE:         score {fscore.score(badD)}, {len(FD)} failures, {len(badD)} bad atoms')
print(f'  failures   : {FD}')

s0, sD = set(F0), set(FD)
print(f'\nfixed by the deliverable (in baseline, not in deliverable): {len(s0-sD)}')
print(f'   {sorted(s0-sD)}')
print(f'NEWLY BROKEN by the deliverable (not in baseline, in deliverable): {len(sD-s0)}')
print(f'   {sorted(sD-s0)}')
print(f'common to both: {len(s0 & sD)} -> {sorted(s0 & sD)}')

print('\n--- the decisive question ---')
print('The greedy tuner only builds rows for equations FAILING AT BASELINE.')
print(f'Of the deliverable\'s {len(FD)} failures, how many are already failing at baseline?'
      f'  {len(sD & s0)}')
print(f'Of the {len(s0-sD)} equations the deliverable FIXES, all are baseline failures by')
print('construction, so they ARE in the row set. So the rows exist -- the issue is that')
print('zeroing them requires ACCEPTING the newly-broken ones, which greedy never models.')

# how big is the delta the deliverable applies?
print('\ndeliverable delta on the freed variables (vs uncorrupted baseline):')
for u in freed:
    d = vd[u] - P.v_unc[u]
    print(f'  x_{u:<6d} delta {len(str(abs(d))) if d else 0} digits')

# is the deliverable's own delta expressible in the affine model built at baseline?
print('\nchecking affinity of the freed handles AT THE BASELINE:')
aff, cols = PR._affine_cols(eng, v0, bad0, sorted(freed))
print(f'  affine among the 5 freed: {sorted(aff)}')
for f in freed:
    if f in cols:
        print(f'  x_{f:<6d} column touches atoms {sorted(cols[f])[:8]}')
    else:
        print(f'  x_{f:<6d} NOT AFFINE at baseline  <-- cannot be modelled linearly here')
