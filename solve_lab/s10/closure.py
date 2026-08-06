"""S10 step 36: the FULL CLOSURE linear system over GF(p).

Iterate to a fixed point:
   columns C  = non-boolean free inputs available as knobs
   rows    R  = every check reachable from C (failing -> -residual, satisfied -> 0)
   expand  C := C + gradient supports of all rows
until stable, then solve over GF(p).

Also diagnose degenerate (square) checks: a check that is 0 with an identically
zero gradient is locally flat, so its linear row is vacuous and the real
constraint is on its square root.  Report any such row rather than silently
trusting it (this was Session 9's methodological bug).
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
MAXCOL = int(os.environ.get('MAXCOL', 3000))

v = L.load(os.path.join(HERE, 'forward_state.json'))
vm = [x % P for x in v]
av = L.all_atom_values(v)
FAIL = set(a for a in range(L.NA) if av[a] and a not in atom_out)
print('failing checks:', sorted(FAIL), flush=True)

deg = [max((len(m) for m in Pp), default=0) for Pp in L.polys]
gcache = {}


def gradc(c):
    if c not in gcache:
        gcache[c] = {u: d for u, d in ad.grad(c, vm).items() if u not in BOOL}
    return gcache[c]


def live_reach_checks(u):
    seen, frontier, checks = {u}, [u], set()
    while frontier:
        nxt = []
        for w in frontier:
            for a in L.var_atoms[w]:
                if ad.dpart(a, w, vm) == 0:
                    continue
                if a not in atom_out:
                    checks.add(a); continue
                t = atom_out[a][1]
                if t == w or t in seen or ad.dpart(a, t, vm) == 0:
                    continue
                seen.add(t); nxt.append(t)
        frontier = nxt
    return checks


cols = set()
for c in FAIL:
    cols |= set(gradc(c))
reach_cache = {}
t0 = time.time()
for rnd in range(8):
    rows_set = set()
    for u in cols:
        if u not in reach_cache:
            reach_cache[u] = live_reach_checks(u)
        rows_set |= reach_cache[u]
    newcols = set(cols)
    for c in rows_set:
        newcols |= set(gradc(c))
    print(f'round {rnd}: cols={len(cols)} rows={len(rows_set)} '
          f'-> newcols={len(newcols)}  ({time.time()-t0:.0f}s)', flush=True)
    if len(newcols) == len(cols) or len(newcols) > MAXCOL:
        cols = newcols
        break
    cols = newcols

rows_set = set()
for u in cols:
    if u not in reach_cache:
        reach_cache[u] = live_reach_checks(u)
    rows_set |= reach_cache[u]
COLS = sorted(cols)
ROWS = sorted(rows_set)
print(f'\nFINAL system: {len(ROWS)} rows x {len(COLS)} columns', flush=True)

# degeneracy diagnostic
flat = []
for c in ROWS:
    if av[c] == 0 and not gradc(c) and deg[c] >= 3:
        flat.append(c)
print(f'degenerate (flat, degree>=3, value 0) rows: {len(flat)} {flat[:10]}', flush=True)

rows, rhs, names = [], [], []
for c in ROWS:
    g = {u: d for u, d in gradc(c).items() if u in cols}
    if not g and (av[c] % P) == 0:
        continue
    rows.append(g); rhs.append((-av[c]) % P); names.append(c)
print(f'non-trivial rows: {len(rows)}', flush=True)

sol = solve_gfp(rows, rhs, COLS)
if sol is None:
    print('\n*** CLOSURE SYSTEM INCONSISTENT', flush=True)
    keep, krhs, kn = [], [], []
    for r, b, nm in zip(rows, rhs, names):
        keep.append(r); krhs.append(b); kn.append(nm)
        if solve_gfp(keep, krhs, COLS) is None:
            print(f'  first inconsistency on adding check a{nm}')
            print(f'  minimal prefix: {kn}')
            break
else:
    print(f'\n*** CONSISTENT -- step over {len(sol)} inputs', flush=True)
    for u, d in sol.items():
        v[u] = v[u] + d
    ad.fwd(v)
    av2 = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av2[a]]
    fail = L.failing_eqs(av2)
    print(f'   after step: nz={nz}\n   failing={len(fail)} score={L.NEQ-len(fail)}')
    T.save(v, os.path.join(HERE, 'closure_step.json'))
