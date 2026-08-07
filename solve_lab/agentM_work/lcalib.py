"""The calibration ladder for the lattice-column pricer.

Same discipline as calib_r14.py, and deliberately the same first rung: an instrument nobody can
trust is not worth building.  The new rungs are the ones that test the NEW code rather than the
inherited engine -- L4 checks that the lattice module's exact columns agree integer-for-integer
with `ieng.affine_cols` composed with the equation coefficient maps, and L5 checks that the new
selector, restricted to the OLD knob set, reproduces the old answer.
"""
import sys, os, json, time, collections

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
os.chdir(MDIR)
sys.path.insert(0, MDIR)
import shim                                                    # noqa: F401
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H                                            # noqa: E402
import ieng, fscore                                            # noqa: E402
import lattice as LT                                           # noqa: E402

D4 = LT.D4
COF12 = [105, 1329, 3387, 5081, 5676, 9413, 10903, 11436, 14393, 14768, 17325, 22820]
T_LIST = [2554, 6816, 8124, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125]
VD = ieng.VD
ok = {}

print(f'frame U baseline : {fscore.score(ieng.BAD_UNC)}, {len(ieng.FAILS_UNC)} failing', flush=True)
print(f'frame D (witness): {LT.SCORE_D}, {len(LT.FAILS_D)} failing  {LT.FAILS_D}', flush=True)

# ---------------- L1: the deliverable, bit for bit, in the new module's own frame ----------
print('\n=== L1: frame D IS the deliverable, byte for byte ===', flush=True)
nd = sum(1 for i in range(LT.NV) if LT.V_D[i] != VD[i])
ok['L1'] = (LT.SCORE_D == 39026 and LT.FAILS_D == LT.DELIV_FAILS and nd == 0
            and len(LT.BAD_D) == 8)
print(f'  score {LT.SCORE_D}  fails {LT.FAILS_D}  nonzero atoms {len(LT.BAD_D)}  '
      f'vars differing {nd} of {LT.NV}', flush=True)
print('  L1', 'PASSED' if ok['L1'] else 'FAILED', flush=True)

# ---------------- L2: the three CLI-agreeing points -----------------------------------
print('\n=== L2: the three CLI-agreeing points ===', flush=True)
g = True
for leaf, exp in ((None, 39026), (4287, 39000), (17378, 38961)):
    ch = {u: VD[u] for u in LT.FREED}
    if leaf is not None:
        ch[leaf] = 1
    sc, _b, _v = ieng.score_from_unc(ch, LT.PIN)
    g &= (sc == exp)
    print(f'  {"deliverable" if leaf is None else f"+x_{leaf}=1":16s} -> {sc}  expect {exp}',
          flush=True)
ok['L2'] = g
print('  L2', 'PASSED' if g else 'FAILED', flush=True)

# ---------------- L3: T's cofactor calibration ----------------------------------------
print('\n=== L3: T calibration (12 cofactors zeroed) ===', flush=True)
ch = {u: VD[u] for u in LT.FREED}
for c in COF12:
    ch[c] = 0
sc, bad, _v = ieng.score_from_unc(ch, LT.PIN)
fl = sorted(fscore.fails(bad))
ok['L3'] = (sc == 39021 and fl == T_LIST)
print(f'  score {sc} (expect 39021)  failing {len(fl)} (expect 12)  list match {fl == T_LIST}',
      flush=True)
print('  L3', 'PASSED' if ok['L3'] else 'FAILED', flush=True)

# ---------------- L4: the NEW columns == the OLD columns, integer for integer -----------
print('\n=== L4: lattice columns == ieng.affine_cols o CM, exactly ===', flush=True)
aff, acols = ieng.affine_cols(LT.PIN, LT.FREED)
CMU = LT.eqmaps(ieng.FAILS_UNC)
mismatch = 0; checked = 0
for u in aff:
    old = {}
    for e in ieng.FAILS_UNC:
        cm, _c0 = CMU[e]
        s = 0
        for a, d in acols[u].items():
            c = cm.get(a)
            if c:
                s += c * d
        if s:
            old[e] = s
    new = LT.column(u, ieng.V_UNC, ieng.BAD_UNC, ieng.FAILS_UNC, CMU, LT.PIN)
    checked += 1
    if new != old:
        mismatch += 1
        print(f'  MISMATCH x_{u}', flush=True)
ok['L4'] = (mismatch == 0 and checked == len(aff))
print(f'  {checked} knob columns compared over {len(ieng.FAILS_UNC)} rows, mismatches {mismatch}',
      flush=True)
print('  L4', 'PASSED' if ok['L4'] else 'FAILED', flush=True)

# ---------------- L5: the new solver on the OLD knob set reproduces the old answer -------
print('\n=== L5: lattice.price on the old knob set == ieng.tune ===', flush=True)
colsU = {u: LT.column(u, ieng.V_UNC, ieng.BAD_UNC, ieng.FAILS_UNC, CMU, LT.PIN) for u in aff}
colsU = {u: c for u, c in colsU.items() if c is not None}
rhsU = {}
for e in ieng.FAILS_UNC:
    cm, c0 = CMU[e]
    rhsU[e] = -(c0 + sum(c * ieng.BAD_UNC[a] for a, c in cm.items() if a in ieng.BAD_UNC))
pr = LT.price(sorted(colsU), colsU, rhsU, ieng.FAILS_UNC, ieng.V_UNC, ieng.BAD_UNC, LT.PIN)
tu = ieng.tune(D4, nprobe=80, budget=180.0)
ok['L5'] = (pr['score'] == tu['score'] == 39026)
print(f'  lattice.price {pr["score"]}   ieng.tune {tu["score"]}   expect 39026 both', flush=True)
print('  L5', 'PASSED' if ok['L5'] else 'FAILED', flush=True)

# ---------------- L6: the lattice diagnostics reproduce a known fact --------------------
print('\n=== L6: SNF diagnostics, frame U, the old knob set ===', flush=True)
rep = LT.snf_report(colsU, sorted(colsU), rhsU, ieng.FAILS_UNC)
print(f'  {rep}', flush=True)
ok['L6'] = (rep['rank_A'] > 0)
print('  L6', 'PASSED' if ok['L6'] else 'FAILED', flush=True)

print('\nLADDER:', {k: ('PASS' if v else 'FAIL') for k, v in ok.items()}, flush=True)
print('ALL PASSED' if all(ok.values()) else 'SOME FAILED', flush=True)
json.dump({k: bool(v) for k, v in ok.items()}, open('lcalib.json', 'w'), indent=1)
