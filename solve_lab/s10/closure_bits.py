"""S10 step 38: same closure, but with the BOOLEAN gate-bits admitted as GF(p)
unknowns.

If the system becomes CONSISTENT once bits are relaxed to GF(p), the obstruction
is not algebraic -- it is the 0/1 integrality of the message bits, i.e. a
combinatorial problem, and the relaxed solution tells us which bits matter.
If it stays inconsistent, the obstruction is algebraic and no bit pattern helps.
"""
import os, sys, collections, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad
from newton import BOOL, solve_gfp

P = ad.P
atom_out = L.atom_out
MAXCOL = int(os.environ.get('MAXCOL', 4000))
USE_BOOL = os.environ.get('NOBOOL', '') == ''

v = L.load(os.path.join(HERE, 'forward_state.json'))
vm = [x % P for x in v]
av = L.all_atom_values(v)
FAIL = set(a for a in range(L.NA) if av[a] and a not in atom_out)
print(f'failing checks: {sorted(FAIL)}   (bits {"INCLUDED" if USE_BOOL else "excluded"})',
      flush=True)

gcache = {}
def gradc(c):
    if c not in gcache:
        g = ad.grad(c, vm)
        if not USE_BOOL:
            g = {u: d for u, d in g.items() if u not in BOOL}
        gcache[c] = g
    return gcache[c]

reach_cache = {}
def live_reach_checks(u):
    if u in reach_cache: return reach_cache[u]
    seen, frontier, checks = {u}, [u], set()
    while frontier:
        nxt = []
        for w in frontier:
            for a in L.var_atoms[w]:
                if ad.dpart(a, w, vm) == 0: continue
                if a not in atom_out:
                    checks.add(a); continue
                t = atom_out[a][1]
                if t == w or t in seen or ad.dpart(a, t, vm) == 0: continue
                seen.add(t); nxt.append(t)
        frontier = nxt
    reach_cache[u] = checks
    return checks

cols = set()
for c in FAIL: cols |= set(gradc(c))
t0 = time.time()
for rnd in range(10):
    rows_set = set()
    for u in cols: rows_set |= live_reach_checks(u)
    newcols = set(cols)
    for c in rows_set: newcols |= set(gradc(c))
    print(f'round {rnd}: cols={len(cols)} rows={len(rows_set)} -> {len(newcols)} '
          f'({time.time()-t0:.0f}s)', flush=True)
    if len(newcols) == len(cols) or len(newcols) > MAXCOL:
        cols = newcols; break
    cols = newcols

rows_set = set()
for u in cols: rows_set |= live_reach_checks(u)
COLS = sorted(cols); ROWS = sorted(rows_set)
nb = len([u for u in COLS if u in BOOL])
print(f'\nsystem: {len(ROWS)} rows x {len(COLS)} cols ({nb} boolean, '
      f'{len(COLS)-nb} non-boolean)', flush=True)

rows, rhs, names = [], [], []
for c in ROWS:
    g = {u: d for u, d in gradc(c).items() if u in cols}
    if not g and (av[c] % P) == 0: continue
    rows.append(g); rhs.append((-av[c]) % P); names.append(c)
print(f'non-trivial rows: {len(rows)}', flush=True)

sol = solve_gfp(rows, rhs, COLS)
print('\n*** RESULT:', 'CONSISTENT' if sol is not None else 'INCONSISTENT', flush=True)
if sol is not None:
    bits = {u: d for u, d in sol.items() if u in BOOL}
    print(f'   solution touches {len(sol)} columns, {len(bits)} of them boolean')
    print(f'   boolean deltas that are NOT 0/+-1: '
          f'{[(f"x_{u}", str(d)[:20]) for u, d in list(bits.items())[:10]]}')
    json.dump({f'x_{u}': str(d) for u, d in sol.items()},
              open(os.path.join(HERE, 'relaxed_sol.json'), 'w'))
else:
    # minimal certificate: left null vector over the failing rows
    print('   searching for a minimal inconsistent subset...')
    import itertools
    idx = list(range(len(rows)))
    cur = idx[:]
    changed = True
    while changed:
        changed = False
        for i in list(cur):
            trial = [j for j in cur if j != i]
            if solve_gfp([rows[j] for j in trial], [rhs[j] for j in trial], COLS) is None:
                cur = trial; changed = True
    print(f'   MINIMAL inconsistent set ({len(cur)} rows): {[names[j] for j in cur]}')
    print(f'   of which failing: {[names[j] for j in cur if names[j] in FAIL]}')
