"""Calibrate the incremental engine before running anything at scale.

Gates, all of which the previous engine passes:
  G1  reproduce the deliverable from its four handles: 39,026 / 7 failing / 8 atoms / 0 vars differing
  G2  the three points my scorer agrees with the checker.py CLI on:
        deliverable                      39,026
        deliverable + x_4287 = 1         39,000   (CLI-verified earlier)
        deliverable + x_17378 = 1        38,961   (CLI-verified earlier)
  G3  T's calibration: deliverable with the 12 cofactors zeroed -> 39,021 / 12 / exact list
  G4  incremental result identical to a full engine3 forward + badatoms (exactness)
  G5  tune() from the shared uncorrupted baseline recovers 39,026 on the deliverable's four
  G6  timing on GENERAL 4-subsets (the case that was slow), vs the old engine
"""
import sys, os, json, time, itertools, pickle
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import ieng, price as PR, fscore
import engine3 as E3

D4 = [642, 28730, 29854, 31864]
COF12 = [105, 1329, 3387, 5081, 5676, 9413, 10903, 11436, 14393, 14768, 17325, 22820]
T_LIST = [2554, 6816, 8124, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125]
D_LIST = [12231, 12270, 12350, 14584, 18673, 22044, 29125]
VD = ieng.VD

print(f'shared baseline: score {fscore.score(ieng.BAD_UNC)}, '
      f'{len(ieng.FAILS_UNC)} failing, {len(ieng.BAD_UNC)} bad atoms', flush=True)

# ---------- G1 ----------
print('\n=== G1: deliverable from its four handles ===', flush=True)
t0 = time.time()
r = ieng.price_given(D4, {u: VD[u] for u in ieng.site(D4)[0]})
dt = time.time() - t0
nd = sum(1 for i in range(ieng.NV) if r['v'][i] != VD[i])
print(f'  freed {r["freed"]}', flush=True)
print(f'  score {r["score"]}  fails {r["fails"]}  nbad {len(r["bad"])}  '
      f'vars differing {nd}   ({dt:.2f}s)', flush=True)
g1 = (r['score'] == 39026 and r['fails'] == D_LIST and nd == 0 and len(r['bad']) == 8)
print(f'  G1 {"PASSED" if g1 else "FAILED"}', flush=True)

# ---------- G2 ----------
print('\n=== G2: the three CLI-agreeing points ===', flush=True)
freed, pin = ieng.site(D4)
base_ch = {u: VD[u] for u in freed}
g2 = True
for extra, expect in ((None, 39026), (4287, 39000), (17378, 38961)):
    ch = dict(base_ch)
    if extra is not None:
        ch[extra] = 1
    sc, bad, v = ieng.score_from_unc(ch, pin)
    ok = (sc == expect)
    g2 &= ok
    lab = 'deliverable' if extra is None else f'deliverable + x_{extra}=1'
    print(f'  {lab:28s} -> {sc}  expect {expect}  {"OK" if ok else "MISMATCH"}', flush=True)
print(f'  G2 {"PASSED" if g2 else "FAILED"}', flush=True)

# ---------- G3 ----------
print("\n=== G3: T's calibration (12 cofactors zeroed) ===", flush=True)
ch = dict(base_ch)
for u in COF12:
    ch[u] = 0
sc3, bad3, v3 = ieng.score_from_unc(ch, pin)
f3 = sorted(fscore.fails(bad3))
g3 = (sc3 == 39021 and f3 == T_LIST)
print(f'  score {sc3} (expect 39021)  failing {len(f3)} (expect 12)  list match {f3 == T_LIST}',
      flush=True)
print(f'  G3 {"PASSED" if g3 else "FAILED"}', flush=True)

# ---------- G4 ----------
print('\n=== G4: incremental == full engine3 forward ===', flush=True)
_, demote = PR.closure(D4)
eng = E3.Eng(demote)
seed = {f: VD[f] for f in eng.FREE if VD[f] != 0}
seed[4287] = 1
vfull = eng.forward(seed)
badfull = eng.badatoms(vfull)
scfull = fscore.score(badfull)
ch = dict(base_ch); ch[4287] = 1
scinc, badinc, vinc = ieng.score_from_unc(ch, pin)
same_v = sum(1 for i in range(ieng.NV) if vfull[i] != vinc[i])
g4 = (scfull == scinc and same_v == 0 and sorted(badfull) == sorted(badinc))
print(f'  full {scfull}  incremental {scinc}  vars differing {same_v}  '
      f'atoms identical {sorted(badfull)==sorted(badinc)}', flush=True)
print(f'  G4 {"PASSED" if g4 else "FAILED"}', flush=True)

# ---------- G5 ----------
print('\n=== G5: tune() recovers 39,026 from the shared baseline ===', flush=True)
t0 = time.time()
rt = ieng.tune(D4, want=True)
print(f'  base {rt["base_score"]} -> {rt["score"]}  ({rt["secs"]:.2f}s, '
      f'{rt["nrows_target"]} rows, {rt["nknobs"]} knobs)', flush=True)
g5 = rt['score'] >= 39026
print(f'  G5 {"PASSED" if g5 else "FAILED"}', flush=True)

# ---------- G6 ----------
print('\n=== G6: timing on GENERAL 4-subsets (the slow case) ===', flush=True)
prev = pickle.load(open('pricelead.pkl', 'rb'))
rtm = {v['added']: v['rt'] for v in prev.values()}
POOL = sorted(rtm, key=lambda u: (-rtm[u], u))
sites = list(itertools.combinations(POOL[:14], 4))[:25]
t0 = time.time(); n = 0
for s in sites:
    ieng.tune(list(s))
    n += 1
per = (time.time() - t0) / n
print(f'  {n} general 4-subsets, mean {per:.3f}s/site', flush=True)
print(f'  old engine measured ~0.4s/site only for SHARED-demotion sites, and 1-2 orders '
      f'slower for general ones', flush=True)
tot = 4251930 * per
print(f'  projected C(102,4)=4,251,930 sites: {tot/3600:.1f} core-hours, '
      f'{tot/3600/4:.1f} h on 4 cores', flush=True)

print(f'\nALL GATES {"PASSED" if (g1 and g2 and g3 and g4 and g5) else "FAILED"}', flush=True)
