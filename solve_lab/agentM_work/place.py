"""Price ALTERNATIVE DEFECT PLACEMENTS by what actually cancels.

The eqsub table priced *repairs of the current placement*. This prices *different
placements*: demote an extra atom (score-neutral, it just frees a variable), then use
that new freedom to satisfy failing equations BY CANCELLATION -- two nonzero atoms in
one equation whose contributions sum to zero -- rather than by driving atoms to zero.

That is the move the deliverable itself makes: it holds 8 atoms nonzero and is seen by
only 7 equations, whereas E's 2-nonzero-atom state fails 28. More nonzero atoms, fewer
failures, because they cancel.

Every candidate is scored the same way as eqsub: solve, APPLY, re-propagate, measure.
"""
import sys, os, json, time, math, itertools, collections, pickle
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as E_broken
import engine3 as E3, mcore2 as M, fscore, sparse

NEQ = len(H.eqt)
vd = M.load_vec()
BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 2100
MAXEXTRA = int(sys.argv[2]) if len(sys.argv) > 2 else 1

cd = pickle.load(open('placecands.pkl', 'rb'))
CANDS = cd['cands']
base_eng = E3.Eng(E3.BASE_DEMOTE)
v_base = base_eng.forward(base_eng.seed_of(vd))
bad_base = base_eng.badatoms(v_base)
FAILS = sorted(fscore.fails(bad_base))
print(f'baseline score {fscore.score(bad_base)}, failing {FAILS}', flush=True)
print(f'{len(CANDS)} demotable candidates: {CANDS}', flush=True)


def resid_delta(eng, v0, base_bad, changes):
    v, aff = eng.apply_delta(v0, changes)
    touched = set()
    for u in set(aff) | set(changes):
        touched.update(E3._ATOM_OF[u])
    ns = {'v': v, '__builtins__': {}}
    bad = dict(base_bad)
    for i in touched:
        r = eval(H.acodes[i], ns)
        if r:
            bad[i] = r
        else:
            bad.pop(i, None)
    return bad, v


def affine_cols(eng, v0, bad0, cand):
    cols = {}; aff = []
    for f in cand:
        o = v0[f]
        try:
            b1, _ = resid_delta(eng, v0, bad0, {f: o + 1})
            b2, _ = resid_delta(eng, v0, bad0, {f: o + 2})
            b7, _ = resid_delta(eng, v0, bad0, {f: o + 7})
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
    return aff, cols


CM = {}
for e in FAILS:
    cm = collections.defaultdict(int); const = 0
    for c, a in H.eqt[e][2]:
        if a < 0:
            const += c
        else:
            cm[a] += c
    CM[e] = (dict(cm), const)

results = {}
best = (fscore.score(bad_base), None, None)
t0 = time.time()

extras = [()] + [(a,) for a in CANDS]
if MAXEXTRA >= 2:
    extras += list(itertools.combinations(CANDS, 2))

for extra in extras:
    if time.time() - t0 > BUDGET:
        print('  [budget]', flush=True); break
    D = list(E3.BASE_DEMOTE) + list(extra)
    eng = E3.Eng(D)
    seed = eng.seed_of(vd)
    v0 = eng.forward(seed)
    nd = sum(1 for u in range(E3.NV) if v0[u] != vd[u])
    bad0 = eng.badatoms(v0)
    sc0 = fscore.score(bad0)
    if nd != 0 or sc0 != 39026:
        print(f'extra={extra}: NOT score-neutral (diff {nd}, score {sc0}) -- skipped', flush=True)
        continue
    newvars = [E3.ATOM2VAR[a] for a in extra]
    # knob set: freed vars (incl. the new ones) + cone of every atom of the failing eqs
    cand = set(eng.pin)
    for e in FAILS:
        for c, a in H.eqt[e][2]:
            if a >= 0:
                try:
                    cand |= set(E_broken.cone(a)[1])
                except Exception:
                    pass
    FS = set(eng.FREE)
    cand = sorted(f for f in cand if f in FS)
    aff, cols = affine_cols(eng, v0, bad0, cand)

    def row_for(e):
        cm, const = CM[e]
        row = {}
        for f in aff:
            co = 0
            for a, d in cols[f].items():
                c = cm.get(a)
                if c:
                    co += c * d
            if co:
                row[f] = co
        s0 = const + sum(c * bad0[a] for a, c in cm.items() if a in bad0)
        return row, -s0

    ROWS = {e: row_for(e) for e in FAILS}
    localbest = (sc0, None)
    for k in range(1, len(FAILS) + 1):
        for S in itertools.combinations(FAILS, k):
            rows = [ROWS[e][0] for e in S]; rhs = [ROWS[e][1] for e in S]
            if any(not r for r in rows):
                continue
            sol, msg, _ = sparse.solve_sparse(rows, rhs, names=list(S), verbose=False,
                                              maxcore=400, maxcorebits=5_000_000)
            if sol is None:
                continue
            ns = dict(seed)
            for f, d in sol.items():
                if d:
                    ns[f] = v0[f] + d
            try:
                v = eng.forward(ns)
                av = eng.badatoms(v)
                sc = fscore.score(av)
            except Exception:
                continue
            if sc > localbest[0]:
                localbest = (sc, S)
            if sc > best[0]:
                best = (sc, ns, (extra, S))
                json.dump({f"x_{i}": int(v[i]) for i in range(E3.NV) if v[i] != 0},
                          open(f'M_place_{sc}.json', 'w'))
                print(f'  *** ABOVE BASELINE {sc}  extra={extra} S={S} -> M_place_{sc}.json ***',
                      flush=True)
        if time.time() - t0 > BUDGET:
            break
    results[extra] = localbest
    print(f'extra={extra} (new knobs {newvars}, {len(aff)} affine): best {localbest[0]} '
          f'via {localbest[1]}   [{time.time()-t0:.0f}s]', flush=True)

pickle.dump(results, open('place.pkl', 'wb'))
print('\n=== SUMMARY ===')
for extra, (sc, S) in sorted(results.items(), key=lambda kv: -kv[1][0]):
    print(f'  extra={extra}: {sc}  via {S}')
print(f'\nBEST OVERALL {best[0]}   baseline 39026')
