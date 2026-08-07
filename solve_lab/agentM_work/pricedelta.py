"""Price O's lattice target: make the whole defect region hold, then measure the collateral.

O inverted rather than sampled and found the support {a23616, a23618, a36660, a36662}
solves, making all 13 region equations hold.  Two carriers are free and zero-collateral
(x_8731 for a36662, x_9118 for a36660); the two open ones are the defining expressions of
handles x642 and x28730 -- two of the four in my calibration row.

I do not have O's integer shift vector, so I re-derive the target in my own frame: solve
the region equations exactly with a knob set that INCLUDES the free carriers.  That is the
one thing my earlier eqsub run lacked -- it used only the freed handles and reported 38,989
for "fix all 7".  Then APPLY and re-propagate: the collateral is measured, not modelled.
"""
import sys, os, json, time, collections, itertools, pickle
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as EB, engine3 as E3
import price as PR, fscore, sparse, sweep

NEQ = PR.NEQ
D4 = [642, 28730, 29854, 31864]
CARRIERS = [8731, 9118, 4432]
REGION13 = json.load(open('lcrit.json'))['L_style_baseline_fails']
FAIL7 = [12231, 12270, 12350, 14584, 18673, 22044, 29125]
OATOMS = [23616, 23618, 36660, 36662]

vd = PR.load_deliverable()
freed, demote = PR.closure(D4)
eng = E3.Eng(demote)
pinset = set(eng.pin)
seed = {f: vd[f] for f in eng.FREE if vd[f] != 0}
v0 = eng.forward(seed)
bad0 = eng.badatoms(v0)
print(f'START = the deliverable: score {fscore.score(bad0)}, bad {sorted(bad0)}', flush=True)
print(f'freed {freed}', flush=True)

FS = set(eng.FREE)
cand = set(freed) | {c for c in CARRIERS if c in FS}
for e in REGION13:
    for c, a in H.eqt[e][2]:
        if a >= 0:
            try:
                cand |= set(EB.cone(a)[1])
            except Exception:
                pass
for a in OATOMS:
    try:
        cand |= set(EB.cone(a)[1])
    except Exception:
        pass
cand = sorted(f for f in cand if f in FS)
print(f'{len(cand)} candidate knobs (carriers included: '
      f'{ {c: (c in cand) for c in CARRIERS} })', flush=True)

t0 = time.time()
cols = {}; aff = []
for f in cand:
    o = v0[f]
    try:
        b1, _ = sweep.fast_resid(eng, v0, bad0, {f: o + 1}, pinset)
        b2, _ = sweep.fast_resid(eng, v0, bad0, {f: o + 2}, pinset)
        b7, _ = sweep.fast_resid(eng, v0, bad0, {f: o + 7}, pinset)
    except Exception:
        continue
    col = {}; ok = True
    for a in set(b1) | set(b2) | set(b7) | set(bad0):
        d1 = b1.get(a, 0) - bad0.get(a, 0)
        if b2.get(a, 0) - bad0.get(a, 0) != 2 * d1 or b7.get(a, 0) - bad0.get(a, 0) != 7 * d1:
            ok = False; break
        if d1:
            col[a] = d1
    if ok:
        aff.append(f); cols[f] = col
print(f'{len(aff)} affine knobs ({time.time()-t0:.0f}s); carriers affine: '
      f'{ {c: (c in cols) for c in CARRIERS} }', flush=True)
for c in CARRIERS:
    if c in cols:
        print(f'   x_{c} column -> atoms {sorted(cols[c])}', flush=True)


def rowfor(e):
    cm = collections.defaultdict(int); const = 0
    for c, a in H.eqt[e][2]:
        if a < 0:
            const += c
        else:
            cm[a] += c
    row = {}
    for f in aff:
        co = 0
        for a, d in cols[f].items():
            c = cm.get(a)
            if c:
                co += c * d
        if co:
            row[f] = co
    return row, -(const + sum(c * bad0[a] for a, c in cm.items() if a in bad0))


def try_target(eqs, label):
    rows = []; rhs = []
    for e in eqs:
        r, b = rowfor(e)
        rows.append(r); rhs.append(b)
    idx = [i for i in range(len(rows)) if rows[i] or rhs[i] != 0]
    sol, msg, _ = sparse.solve_sparse([rows[i] for i in idx], [rhs[i] for i in idx],
                                      names=[eqs[i] for i in idx], verbose=False,
                                      maxcore=600, maxcorebits=8_000_000)
    if sol is None:
        print(f'  [{label}] simultaneous solve: NO SOLUTION ({msg})', flush=True)
        return None
    ch = {}
    for f, d in sol.items():
        if d:
            ch[f] = v0[f] + d
    if not ch:
        print(f'  [{label}] solution is trivial (no change)', flush=True)
        return None
    ns = dict(seed); ns.update(ch)
    v = eng.forward(ns)
    av = eng.badatoms(v)
    fl = sorted(fscore.fails(av))
    sc = NEQ - len(fl)
    got = [e for e in eqs if e not in fl]
    print(f'  [{label}] SOLVED. {len(got)}/{len(eqs)} targets now hold. '
          f'score {sc}  (nbad {len(av)}, nfail {len(fl)})', flush=True)
    print(f'      COLLATERAL: {len(fl)} equations still fail -> {fl[:24]}', flush=True)
    print(f'      knobs moved: {len(ch)}  shift sizes (bits): '
          f'{ {f: (v0[f]-ch[f]).bit_length() for f in list(ch)[:8]} }', flush=True)
    if sc > 39026:
        fn = f'M_delta_{sc}.json'
        json.dump({f"x_{k}": int(v[k]) for k in range(PR.NV) if v[k] != 0}, open(fn, 'w'))
        print(f'      *** ABOVE BASELINE -> {fn} ***', flush=True)
    return sc, fl, ns, v


print('\n=== targets ===', flush=True)
res = {}
res['fail7'] = try_target(FAIL7, 'all 7 currently-failing')
res['region13'] = try_target(REGION13, 'all 13 region equations')

# also: every subset of the 7, now WITH the carriers in the knob set
print('\n=== subsets of the 7, with carriers in the knob set ===', flush=True)
best = (39026, None)
for k in range(1, 8):
    for S in itertools.combinations(FAIL7, k):
        rows = []; rhs = []
        for e in S:
            r, b = rowfor(e)
            rows.append(r); rhs.append(b)
        if any(not r for r in rows):
            continue
        sol, _, _ = sparse.solve_sparse(rows, rhs, verbose=False,
                                        maxcore=600, maxcorebits=8_000_000)
        if sol is None:
            continue
        ch = {}
        for f, d in sol.items():
            if d:
                ch[f] = v0[f] + d
        if not ch:
            continue
        try:
            bad, _ = sweep.fast_resid(eng, v0, bad0, ch, pinset)
            sc = fscore.score(bad)
        except Exception:
            continue
        if sc > best[0]:
            best = (sc, S)
            ns = dict(seed); ns.update(ch)
            v = eng.forward(ns)
            fn = f'M_delta_sub_{sc}.json'
            json.dump({f"x_{i}": int(v[i]) for i in range(PR.NV) if v[i] != 0}, open(fn, 'w'))
            print(f'  *** {sc} at subset {S} -> {fn} ***', flush=True)
print(f'best over subsets: {best}', flush=True)
pickle.dump({'best_subset': best}, open('pricedelta.pkl', 'wb'))
print('\nbaseline to beat: 39026 ; perfect: 39033', flush=True)
