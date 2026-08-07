"""Post-restart recalibration (round 14).  Same gates as icalib.py G1-G5, but imports
`shim` FIRST so `harness` resolves to harness_m.py / the model rebuilt in agentM_work.
G6 (timing against pricelead.pkl) is dropped: that pkl was wiped and it was a timing
reference, not a correctness gate.
"""
import sys, os, json, time, itertools

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
os.chdir(MDIR)
sys.path.insert(0, MDIR)
import shim                                                    # noqa: F401
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H                                            # noqa: E402
import ieng, price as PR, fscore                               # noqa: E402
import engine3 as E3                                           # noqa: E402

D4 = [642, 28730, 29854, 31864]
COF12 = [105, 1329, 3387, 5081, 5676, 9413, 10903, 11436, 14393, 14768, 17325, 22820]
T_LIST = [2554, 6816, 8124, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125]
D_LIST = [12231, 12270, 12350, 14584, 18673, 22044, 29125]
VD = ieng.VD

print(f'model: atoms {len(H.atoms)}  eqs {len(H.eqt)}  free {len(H.FREE)}  seq {len(H.SEQ)}',
      flush=True)
print(f'shared baseline: score {fscore.score(ieng.BAD_UNC)}, '
      f'{len(ieng.FAILS_UNC)} failing, {len(ieng.BAD_UNC)} bad atoms', flush=True)

ok = {}

# ---------- G1 ----------
print('\n=== G1: deliverable from its four handles ===', flush=True)
t0 = time.time()
r = ieng.price_given(D4, {u: VD[u] for u in ieng.site(D4)[0]})
dt = time.time() - t0
nd = sum(1 for i in range(ieng.NV) if r['v'][i] != VD[i])
print(f'  freed {r["freed"]}', flush=True)
print(f'  score {r["score"]}  fails {r["fails"]}  nbad {len(r["bad"])}  '
      f'vars differing {nd} of {ieng.NV}   ({dt:.2f}s)', flush=True)
ok['G1'] = (r['score'] == 39026 and r['fails'] == D_LIST and nd == 0 and len(r['bad']) == 8)
print('  G1', 'PASSED' if ok['G1'] else 'FAILED', flush=True)

# ---------- G2 ----------
print('\n=== G2: the three CLI-agreeing points ===', flush=True)
g2 = True
for leaf, exp in ((None, 39026), (4287, 39000), (17378, 38961)):
    vals = {u: VD[u] for u in ieng.site(D4)[0]}
    ch = dict(vals)
    if leaf is not None:
        ch[leaf] = 1
    freed, pin = ieng.site(D4)
    sc, bad, v = ieng.score_from_unc(ch, pin)
    lbl = 'deliverable' if leaf is None else f'deliverable + x_{leaf}=1'
    good = (sc == exp)
    g2 &= good
    print(f'  {lbl:28s} -> {sc}  expect {exp}  {"OK" if good else "MISMATCH"}', flush=True)
ok['G2'] = g2
print('  G2', 'PASSED' if g2 else 'FAILED', flush=True)

# ---------- G3 ----------
print('\n=== G3: T calibration (12 cofactors zeroed) ===', flush=True)
freed, pin = ieng.site(D4)
ch = {u: VD[u] for u in freed}
for c in COF12:
    ch[c] = 0
sc, bad, v = ieng.score_from_unc(ch, pin)
fl = sorted(fscore.fails(bad))
ok['G3'] = (sc == 39021 and fl == T_LIST)
print(f'  score {sc} (expect 39021)  failing {len(fl)} (expect 12)  list match {fl == T_LIST}',
      flush=True)
print('  G3', 'PASSED' if ok['G3'] else 'FAILED', flush=True)

# ---------- G4: incremental == full engine3 ----------
print('\n=== G4: incremental == full engine3 forward ===', flush=True)
freed, pin = ieng.site(D4)
_fr, demote = PR.closure(D4)
E = E3.Eng(demote)
seed = {f: VD[f] for f in E.FREE if VD[f] != 0}
seed[4287] = 1
vfull = E.forward(seed)
badf = E.badatoms(vfull)
scf = fscore.score(badf)
chi = {u: VD[u] for u in freed}
chi[4287] = 1
sci, badi, vi = ieng.score_from_unc(chi, pin)
nd4 = sum(1 for i in range(ieng.NV) if vi[i] != vfull[i])
ok['G4'] = (scf == sci and nd4 == 0 and set(badf) == set(badi))
print(f'  full {scf}  incremental {sci}  vars differing {nd4}  '
      f'atoms identical {set(badf) == set(badi)}', flush=True)
print('  G4', 'PASSED' if ok['G4'] else 'FAILED', flush=True)

# ---------- G5: tune ----------
print('\n=== G5: tune() recovers 39,026 from the shared baseline ===', flush=True)
t0 = time.time()
tr = ieng.tune(D4)
ok['G5'] = (tr['score'] == 39026)
print(f'  base {tr["base_score"]} -> {tr["score"]}  ({tr["secs"]:.2f}s, '
      f'{tr["nrows_target"]} rows, {tr["nknobs"]} knobs)', flush=True)
print('  G5', 'PASSED' if ok['G5'] else 'FAILED', flush=True)

# ---------- G5b: tune at the raised granularity actually used ----------
print('\n=== G5b: tune() at nprobe=80 budget=180 (the p80 run granularity) ===', flush=True)
tr2 = ieng.tune(D4, nprobe=80, budget=180.0)
ok['G5b'] = (tr2['score'] == 39026)
print(f'  base {tr2["base_score"]} -> {tr2["score"]}  ({tr2["secs"]:.2f}s)', flush=True)
print('  G5b', 'PASSED' if ok['G5b'] else 'FAILED', flush=True)

# ---------- G6: timing on general 4-subsets ----------
print('\n=== G6: timing on GENERAL 4-subsets ===', flush=True)
PF = json.load(open('pfamily.json'))
pool = sorted({v['h'] for v in PF['incident_25'].values()})
t0 = time.time(); n = 0
for c in itertools.islice(itertools.combinations(pool, 4), 40):
    ieng.tune(list(c)); n += 1
print(f'  {n} general 4-subsets in {time.time()-t0:.1f}s  '
      f'= {(time.time()-t0)/n:.3f} s/site', flush=True)

print('\nGATES:', {k: ('PASS' if v else 'FAIL') for k, v in ok.items()}, flush=True)
print('ALL PASSED' if all(ok.values()) else 'SOME FAILED', flush=True)
