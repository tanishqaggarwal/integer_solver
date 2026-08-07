"""CL: the honest system.  Many of the 707 closure columns are BOOLEAN-constrained free
inputs (a check x*x - x pins them to {0,1}); the linear model may not move them.
Restrict to legal columns and price the minimum equation sacrifice."""
import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0 = L.all_atom_values(v0)
BADA = [a for a in range(L.NA) if a not in atom_out and av0[a]]
GAIN = set()
for a in BADA: GAIN |= set(L.atom2eq.get(a, {}))

# ---- boolean-constrained variables: a check that is  c*x*x + d*x  (roots 0 and -d/c)
BOOLV = set()
for a in range(L.NA):
    if a in atom_out: continue
    Pp = L.polys[a]
    vs = set()
    for mo in Pp: vs |= set(mo)
    if len(vs) != 1: continue
    u = next(iter(vs))
    if max(len(mo) for mo in Pp) == 2: BOOLV.add(u)
print(f'variables pinned by a single-variable quadratic check: {len(BOOLV)}')

D = json.load(open(os.path.join(HERE, 'cl_closure_cols.json')))
cols = {int(u): {int(c): int(d) for c, d in mm.items()} for u, mm in D['cols'].items()}
UALL = sorted(cols)
ULEG = [u for u in UALL if u not in BOOLV]
print(f'closure columns {len(UALL)}; boolean-pinned among them {len(UALL)-len(ULEG)}; legal {len(ULEG)}')
atoms = sorted(set().union(*[set(c) for c in cols.values()]))
aset = set(atoms)
EQ = sorted(set().union(*[set(L.atom2eq.get(a, {})) for a in atoms]))

def build(U, drop=()):
    m = len(U)
    rowsE = [e for e in EQ if e not in drop]
    n = len(rowsE)
    A = [[0]*(m+1) for _ in range(n)]
    for i, e in enumerate(rowsE):
        mult, sq, co = L.eq_atoms[e]
        s = 0
        for a, c in co.items():
            s += c*av0[a]
            if a in aset:
                for j, u in enumerate(U):
                    d = cols[u].get(a)
                    if d: A[i][j] = (A[i][j] + c*d) % P
        A[i][m] = (-s) % P
    return A, rowsE

def elim(A, n, m, track=True):
    Tr = [{i: 1} for i in range(n)] if track else None
    r = 0
    for c in range(m):
        pr = next((i for i in range(r, n) if A[i][c]), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]
        if Tr: Tr[r], Tr[pr] = Tr[pr], Tr[r]
        inv = pow(A[r][c], -1, P)
        if inv != 1:
            A[r] = [x*inv % P for x in A[r]]
            if Tr: Tr[r] = {k: x*inv % P for k, x in Tr[r].items()}
        Ar = A[r]
        for i in range(r+1, n):
            f = A[i][c]
            if not f: continue
            A[i] = [(A[i][k]-f*Ar[k]) % P for k in range(m+1)]
            if Tr:
                ti = Tr[i]
                for k, x in Tr[r].items():
                    nv = (ti.get(k, 0)-f*x) % P
                    if nv: ti[k] = nv
                    elif k in ti: del ti[k]
        r += 1
    obs = [i for i in range(r, n) if A[i][m] % P and not any(A[i][j] for j in range(m))]
    return r, obs, Tr

for U, tag in [(UALL, 'ALL 707 columns'), (ULEG, f'LEGAL {len(ULEG)} columns')]:
    t0 = time.time()
    A, rowsE = build(U)
    r, obs, Tr = elim(A, len(rowsE), len(U))
    print(f'\n[{tag}] equation system {len(rowsE)} x {len(U)}: rank {r}, '
          f'inconsistent {len(obs)}  ({time.time()-t0:.0f}s)')
    if obs:
        sups = sorted(({rowsE[k] for k, w in Tr[i].items() if w % P} for i in obs), key=len)
        print(f'   smallest obstruction support: {len(sups[0])} equations {sorted(sups[0])[:14]}')
        inter = set.intersection(*sups) if len(sups) > 1 else sups[0]
        print(f'   equations common to ALL obstructions: {len(inter)} {sorted(inter)[:14]}')
        uni = set().union(*sups)
        print(f'   union of obstruction supports: {len(uni)}')
        print(f'   of the currently-failing {len(GAIN)} eqs, in the union: '
              f'{len(GAIN & uni)}; in the intersection: {len(GAIN & inter)}')
        # greedy: drop equations to kill obstructions
        drop = set()
        for step in range(30):
            A2, rowsE2 = build(U, drop)
            r2, obs2, Tr2 = elim(A2, len(rowsE2), len(U))
            print(f'   greedy step {step}: dropped {len(drop)} eqs '
                  f'({len(drop - GAIN)} beyond the current failures), obstructions {len(obs2)}')
            if not obs2:
                print(f'   *** CONSISTENT after sacrificing {len(drop)} equations '
                      f'-> best reachable score {L.NEQ - len(drop)}')
                json.dump(sorted(drop), open(os.path.join(HERE, f'cl_sacrifice_{len(U)}.json'), 'w'))
                break
            cnt = collections.Counter()
            for i in obs2:
                cnt.update({rowsE2[k] for k, w in Tr2[i].items() if w % P})
            pick = max(cnt, key=lambda e: (e in GAIN, cnt[e]))
            drop.add(pick)
