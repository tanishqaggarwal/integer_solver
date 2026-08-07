"""Frame D: the exact lattice of the deliverable's own residual.

Rows = the deliverable's 7 failing equations.  Knobs = every free input that can reach an atom
of one of them, found by backward reachability rather than by handle incidence.  Columns are
exact integer derivatives verified affine.  Then: what does the lattice say, and what does it
cost to act on it?
"""
import sys, os, json, time, itertools, collections

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
os.chdir(MDIR)
sys.path.insert(0, MDIR)
import shim                                                    # noqa: F401
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H                                            # noqa: E402
import ieng, fscore                                            # noqa: E402
import lattice as LT                                           # noqa: E402
from math import gcd

F = LT.FAILS_D
print(f'frame D: score {LT.SCORE_D}, failing {F}', flush=True)

POOL, ATOMS = LT.reach_pool(F)
print(f'atoms in the 7 equations : {len(ATOMS)}', flush=True)
print(f'REACHABILITY POOL        : {len(POOL)} free inputs  {POOL}', flush=True)
print(f'  old handle knobs inside: {sorted(set(POOL) & set(LT.FREED))}', flush=True)

t0 = time.time()
COLS, DROP = {}, []
for u in POOL:
    c = LT.column(u, LT.V_D, LT.BAD_D, F, LT.CM_D, LT.PIN)
    if c is None:
        DROP.append(u)
    elif c:
        COLS[u] = c
K = sorted(COLS)
print(f'columns computed in {time.time()-t0:.1f}s: {len(K)} nonzero, '
      f'{len(POOL)-len(K)-len(DROP)} zero, {len(DROP)} dropped as NON-AFFINE {DROP}', flush=True)

RHS = {e: -LT.eqval(LT.BAD_D, e) for e in F}
print('\nrow residuals (digits): ' + '  '.join(f'{e}:{len(str(abs(RHS[e])))}' for e in F),
      flush=True)

# ---- which knobs move which row, and is a single row solvable over Z ----
print('\nper-row structure:', flush=True)
single = {}
for e in F:
    movers = [u for u in K if e in COLS[u]]
    g = 0
    for u in movers:
        g = gcd(g, abs(COLS[u][e]))
    solv = bool(g) and RHS[e] % g == 0
    single[e] = {'movers': len(movers), 'gcd_digits': len(str(g)) if g else 0, 'Z_solvable': solv}
    print(f'  eq {e}: {len(movers):2d} movers  gcd {"1" if g == 1 else str(g)[:24] + "..."}  '
          f'single-row Z-solvable {solv}', flush=True)

# ---- the lattice of the whole 7-row system ----
print('\nSNF of the full 7 x %d system:' % len(K), flush=True)
rep = LT.snf_report(COLS, K, RHS, F)
print(f'  {rep}', flush=True)

# ---- and of every subset of rows: which are jointly reachable over Z ----
print('\njoint solvability by row subset (all 127):', flush=True)
byk = collections.Counter(); solvable_sets = []
for r in range(1, len(F) + 1):
    for S in itertools.combinations(F, r):
        if LT._int_solvable(COLS, K, RHS, list(S)):
            byk[r] += 1; solvable_sets.append(S)
for r in range(1, len(F) + 1):
    tot = len(list(itertools.combinations(F, r)))
    print(f'  |S|={r}: {byk[r]:3d} of {tot:3d} Z-solvable', flush=True)

# ---- price the ones that are solvable: solve, APPLY, re-propagate, measure ----
print('\npricing every Z-solvable row subset by re-propagation:', flush=True)
best = (LT.SCORE_D, None)
res = []
t0 = time.time()
for S in solvable_sets:
    pr = LT.price(K, COLS, RHS, list(S), LT.V_D, LT.BAD_D, LT.PIN, want=True)
    res.append({'rows': list(S), 'score': pr['score']})
    if pr['score'] > best[0]:
        best = (pr['score'], S)
        ch = pr.get('changes')
        if ch:
            bad, v = ieng.resid(LT.V_D, LT.BAD_D, ch, LT.PIN)
            fn = f'M_lat_{pr["score"]}_{"-".join(map(str, S))}.json'
            json.dump({f"x_{k}": int(v[k]) for k in range(LT.NV) if v[k] != 0}, open(fn, 'w'))
            print(f'  *** ABOVE 39026: {pr["score"]} buying {S} -> {fn} ***', flush=True)
dist = collections.Counter(r['score'] for r in res)
print(f'  {len(res)} priced in {time.time()-t0:.0f}s', flush=True)
for s in sorted(dist, reverse=True):
    print(f'    {s}: {dist[s]}', flush=True)
print(f'BEST {best[0]}  rows {best[1]}', flush=True)

json.dump({'pool': POOL, 'nonzero_knobs': K, 'dropped_nonaffine': DROP,
           'per_row': {str(k): v for k, v in single.items()},
           'snf': rep, 'solvable_by_size': dict(byk),
           'priced': res, 'best': [best[0], list(best[1]) if best[1] else None]},
          open('lrun.json', 'w'), indent=1)
print('wrote lrun.json', flush=True)
