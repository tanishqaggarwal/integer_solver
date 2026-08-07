"""(A) Re-derive L's exact incidence criterion myself, and (B) explain the baseline gap.

(A) L's criterion: every residual atom has exactly one FREE cofactor u that occurs nowhere
    else, so  equation e contains atom a  <=>  u_a in vars(e),  read off checker's varsets.
    I implement it independently rather than reading L's file, which also checks the claim.

(B) The baseline discrepancy is mine to explain.  My baseline un-corrupts to a fixpoint and
    fails 25.  L's zeroes the 16 tuned handle/cofactor variables and fails 13, a strict
    subset.  Construct L's baseline in my own frame and find out what the extra 12 are.
"""
import sys, os, re, json, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as EB, engine3 as E3
import price as PR, fscore
import checker

MY25 = json.load(open('baseline_sets.json'))['A']
D4 = [642, 28730, 29854, 31864]
# L's 16 tuned variables = the 4 corrupted handles + the 12 that differ mod p
L12 = [105, 1329, 3387, 5081, 5676, 9413, 10903, 11436, 14393, 14768, 17325, 22820]
L16 = sorted(set(D4) | set(L12))

codes, varsets = checker.load_equations()
print(f'checker: {len(varsets)} equations; harness eqt: {len(H.eqt)}  aligned={len(varsets)==len(H.eqt)}')

# ---------- (A) L's criterion, implemented independently ----------
FS = set(EB.FREE)
occlen = {u: len(H.occ.get(u, ())) for u in range(H.NV)}

cofactor = {}       # atom -> its unique free, singly-occurring cofactor
for u in H.SEQ:
    a = H.definer[u][0]
    t = H.atoms[a]
    m = re.fullmatch(r'x_%d - x_(\d+) \* x_(\d+)' % u, t)
    if not m:
        continue
    i, j = int(m.group(1)), int(m.group(2))
    cands = [x for x in (i, j) if x in FS and occlen.get(x, 0) == 1]
    if len(cands) == 1:
        cofactor[a] = (u, cands[0])
print(f'\n(A) residual atoms with exactly one free singly-occurring cofactor: {len(cofactor)}')


def incident_against(targets):
    tv = set()
    for e in targets:
        tv |= set(varsets[e])
    out = {}
    for a, (h, u) in cofactor.items():
        if u in tv:
            rt = sum(1 for e in targets if u in varsets[e])
            out[h] = (a, u, rt)
    return out


inc25 = incident_against(MY25)
print(f'    incident against MY 25-equation baseline : {len(inc25)} handles')
for h, (a, u, rt) in sorted(inc25.items(), key=lambda kv: -kv[1][2]):
    mark = '  <== DELIVERABLE' if h in D4 else ''
    print(f'      x{h:<6d} atom {a:<6d} cofactor x{u:<6d} rt {rt:2d}{mark}')

named = [23754, 35619, 9629, 37413, 34113, 28355]
print(f'\n    cross-check against the handles the coordinator named from L:')
for h in named:
    print(f'      x{h:<6d} in my derivation: {h in inc25}')

pool32 = set(json.load(open('incident_pool.json'))['incident_handles'])
print(f'\n    my earlier 32-pool vs this derivation: '
      f'|derived| {len(inc25)}, in 32 {len(set(inc25) & pool32)}, '
      f'NOT in 32 {sorted(set(inc25) - pool32)}')

# ---------- (B) L's baseline, constructed in my frame ----------
print('\n(B) BASELINE RECONCILIATION')
vd = PR.load_deliverable()
freed, demote = PR.closure(D4)
eng = E3.Eng(demote)

# L's baseline: the deliverable with the 16 tuned variables ZEROED
seedL = {f: vd[f] for f in eng.FREE if vd[f] != 0}
for u in L16:
    seedL.pop(u, None)
vL = eng.forward(seedL)
badL = eng.badatoms(vL)
FL = sorted(fscore.fails(badL))
print(f'  L-style baseline (16 zeroed) : score {fscore.score(badL)}, {len(FL)} failures, '
      f'{len(badL)} bad atoms')
print(f'  my baseline (un-corrupt)     : 39008, {len(MY25)} failures')
sL, sM = set(FL), set(MY25)
print(f'  L-style subset of mine? {sL <= sM}    |L| {len(sL)}  |M| {len(sM)}')
print(f'  in mine not L-style ({len(sM-sL)}): {sorted(sM - sL)}')
print(f'  in L-style not mine ({len(sL-sM)}): {sorted(sL - sM)}')

EXTRA = sorted(sM - sL)
print(f'\n  --- what makes the extra {len(EXTRA)} fail in MY baseline? ---')
vM = EB.forward({f: vd[f] for f in EB.FREE if vd[f] != 0})
badM = EB.badatoms(vM)
print(f'  my baseline bad atoms   : {sorted(badM)}')
print(f'  L-style baseline bad atoms: {sorted(badL)}')
for e in EXTRA[:14]:
    inM = sorted(a for a in badM if a in {x for x in H.avars and []} or False)
    atomsM = sorted(a for c, a in H.eqt[e][2] if a >= 0 and a in badM)
    atomsL = sorted(a for c, a in H.eqt[e][2] if a >= 0 and a in badL)
    print(f'    eq {e:6d}: nonzero atoms in MY baseline {atomsM} | in L-style {atomsL}')

json.dump({'derived_incident': {str(h): {'atom': a, 'cofactor': u, 'rt': rt}
                                for h, (a, u, rt) in inc25.items()},
           'L_style_baseline_fails': FL, 'my_baseline_fails': MY25,
           'extra_in_mine': EXTRA},
          open('lcrit.json', 'w'), indent=1)
print('\nwrote lcrit.json')
