"""Model or solver?  Feed the deliverable's OWN delta into the linear model.

If the linear rows, evaluated at the known-good delta, predict zero for the 18 equations
the deliverable fixes, the model is right and my greedy solver simply is not finding that
solution.  If they do not, the affine model built from +1/+2/+7 probes does not hold at a
delta of ~10^728 and the whole linear tuning approach is invalid at this scale.
"""
import sys, os, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine3 as E3
import price as PR, fscore

vd = PR.load_deliverable()
P = PR.TunedPricer(vd)
h0 = [642, 28730, 29854, 31864]
freed, demote = PR.closure(h0)
eng = E3.Eng(demote)

seed_b = {f: vd[f] for f in eng.FREE if vd[f] != 0}
for u in freed:
    if P.v_unc[u]:
        seed_b[u] = P.v_unc[u]
    else:
        seed_b.pop(u, None)
v0 = eng.forward(seed_b)
bad0 = eng.badatoms(v0)
F0 = sorted(fscore.fails(bad0))
FD = [12231, 12270, 12350, 14584, 18673, 22044, 29125]
FIXED = [e for e in F0 if e not in FD]
print(f'baseline {fscore.score(bad0)}, {len(F0)} failures; deliverable fixes {len(FIXED)}')

# the known-good delta, on the freed vars only
delta = {u: vd[u] - v0[u] for u in freed}
print(f'delta on freed vars: { {u: len(str(abs(d))) for u, d in delta.items()} } (digits)')

# 1) does applying it actually give the deliverable?
ns = dict(seed_b)
for u, d in delta.items():
    ns[u] = v0[u] + d
v = eng.forward(ns)
print(f'applying the delta -> score {fscore.score(eng.badatoms(v))}  '
      f'(vars differing from deliverable {sum(1 for i in range(PR.NV) if v[i] != vd[i])})')

# 2) build the affine model at baseline over just the freed vars, and PREDICT
aff, cols = PR._affine_cols(eng, v0, bad0, sorted(freed))
print(f'affine freed vars: {sorted(aff)}')


def eqval_from_atoms(av, e):
    s = 0
    for c, a in H.eqt[e][2]:
        if a < 0:
            s += c
        elif a in av:
            s += c * av[a]
    return s


agree = dis = 0
for e in FIXED[:12]:
    cm = collections.defaultdict(int); const = 0
    for c, a in H.eqt[e][2]:
        if a < 0:
            const += c
        else:
            cm[a] += c
    s0 = const + sum(c * bad0[a] for a, c in cm.items() if a in bad0)
    pred = s0
    for f in aff:
        d = delta.get(f, 0)
        if not d:
            continue
        co = 0
        for a, cd in cols[f].items():
            c = cm.get(a)
            if c:
                co += c * cd
        pred += co * d
    actual = eqval_from_atoms(eng.badatoms(v), e)
    ok = (pred == actual)
    agree += ok; dis += (not ok)
    print(f'  eq {e:6d}: model predicts {"ZERO" if pred==0 else "nonzero"}, '
          f'actual {"ZERO" if actual==0 else "nonzero"}  match={ok}')
print(f'\nmodel/actual agreement on sampled fixed equations: {agree} agree, {dis} disagree')
print('  -> if agreement is high, the model is sound and the SOLVER is the weak part.')
