"""Is 39,026 optimal in equation space?  Window-free test.

The "obstruction on eq 29125" turned out to be a JOINT effect, not a property of that
row (gcd of its knob coefficients is 1, so it is individually satisfiable).  The window
question then becomes: which SUBSET of the 7 failing equations can be zeroed, and what
does the true score do when you actually do it?

This removes the window choice entirely: for each subset S of the 7 failures we solve
"every row in S becomes 0", APPLY it, re-propagate, and measure the exact score with
fscore.  Collateral damage to every other equation is then counted by measurement, not
by a modelling assumption.
"""
import sys, os, json, time, math, itertools, collections, pickle
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as E_broken
import engine2 as E2, fast2, mcore2 as M, chan2 as C, fscore, sparse

NEQ = len(H.eqt)
vd = M.load_vec()
base = E2.seed_of(vd)
v0 = E2.forward(base)
bad0 = E2.badatoms(v0)
FAILS = sorted(fscore.fails(bad0))
print(f'base score {fscore.score(bad0)}  failing {FAILS}', flush=True)

# ---- widest affine knob set that can move any of the 7 failing equations ----
FREESET = set(E2.FREE)
cand = set(E2.PIN) | set(C.CLUSTERKN)
for e in FAILS:
    for c, a in H.eqt[e][2]:
        if a >= 0:
            try:
                cand |= set(E_broken.cone(a)[1])
            except Exception:
                pass
cand = sorted(f for f in cand if f in FREESET)
aff, cols = C.affine_cols(v0, bad0, cand)
print(f'candidates {len(cand)}, affine knobs {len(aff)}', flush=True)
print(f'the 5 freed definer vars in the knob set: { {f: (f in cols) for f in E2.PIN} }', flush=True)

cmaps = {}
for e in FAILS:
    cm = collections.defaultdict(int)
    const = 0
    for c, a in H.eqt[e][2]:
        if a < 0:
            const += c
        else:
            cm[a] += c
    cmaps[e] = (dict(cm), const)


def row_for(e):
    cm, const = cmaps[e]
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
for e in FAILS:
    r, rhs = ROWS[e]
    g = 0
    for co in r.values():
        g = math.gcd(g, abs(co))
    print(f'  eq {e}: {len(r)} knobs, gcd={g}, individually solvable={bool(g) and (-rhs) % g == 0}',
          flush=True)

results = {}
best = (fscore.score(bad0), None, ())
t0 = time.time()
for k in range(1, len(FAILS) + 1):
    for S in itertools.combinations(FAILS, k):
        rows = [ROWS[e][0] for e in S]
        rhs = [ROWS[e][1] for e in S]
        if any(not r for r in rows):
            results[S] = ('no-knob-dependence', None)
            continue
        sol, msg, _ = sparse.solve_sparse(rows, rhs, names=list(S), verbose=False,
                                          maxcore=400, maxcorebits=5_000_000)
        if sol is None:
            results[S] = ('infeasible', None)
            continue
        ns = dict(base)
        for f, d in sol.items():
            if d:
                ns[f] = v0[f] + d
        try:
            v = E2.forward(ns)
            av = E2.badatoms(v)
            sc = fscore.score(av)
        except Exception as ex:
            results[S] = (f'ERR {type(ex).__name__}', None)
            continue
        nf = sorted(fscore.fails(av))
        got = [e for e in S if e not in nf]
        results[S] = ('ok', sc)
        if sc > best[0]:
            best = (sc, ns, S)
            C.dump(v, f'M_sub_{sc}.json')
            print(f'  *** ABOVE BASELINE {sc} from subset {S} -> M_sub_{sc}.json ***', flush=True)
        print(f'  S={S}: solved, {len(got)}/{len(S)} of them actually zeroed, '
              f'score {sc}  (nbad {len(av)}, nfail {len(nf)})', flush=True)
    print(f'-- size {k} done, {time.time()-t0:.0f}s, best so far {best[0]}', flush=True)

pickle.dump(results, open('eqsub.pkl', 'wb'))
ok = [(S, sc) for S, (st, sc) in results.items() if st == 'ok']
inf = [S for S, (st, sc) in results.items() if st == 'infeasible']
print(f'\n{len(ok)} subsets solvable, {len(inf)} infeasible')
if ok:
    print('best 10 by score:', sorted(ok, key=lambda x: -x[1])[:10])
print('largest solvable subset size:', max((len(S) for S, sc in ok), default=0))
print('BEST OVERALL:', best[0], 'baseline 39026')
