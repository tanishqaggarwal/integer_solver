"""CL: minimum-cost sacrifice set.  Which checks must be released for the cluster's
mod-p system to become consistent, and what does that cost in equations?"""
import os, sys, json, time, collections, itertools
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0 = L.all_atom_values(v0)
BAD = [a for a in range(L.NA) if a not in atom_out and av0[a]]
GAIN = set()
for a in BAD: GAIN |= set(L.atom2eq.get(a, {}))
print(f'failing atoms {BAD}; equations recovered if all vanish: {len(GAIN)}')

D = json.load(open(os.path.join(HERE,'cl_closure_cols.json')))
cols = {int(u): {int(c): int(d) for c, d in m.items()} for u, m in D['cols'].items()}
U = sorted(cols); m = len(U)
allrows = sorted(set().union(*[set(c) for c in cols.values()]) | set(BAD))

def analyse(drop, want_cert=True):
    rows = [a for a in allrows if a not in drop]
    n = len(rows); ri = {a: i for i, a in enumerate(rows)}
    A = [[0]*(m+1) for _ in range(n)]
    for j, u in enumerate(U):
        for a, d in cols[u].items():
            if a in ri: A[ri[a]][j] = d % P
    for a in rows: A[ri[a]][m] = (-av0[a]) % P
    Tr = [{i: 1} for i in range(n)] if want_cert else None
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
    certs = []
    if Tr:
        for i in sorted(obs, key=lambda i: len(Tr[i]))[:14]:
            certs.append({rows[k] for k, w in Tr[i].items() if w % P})
    return len(obs), certs

t0 = time.time()
nobs, certs = analyse(set())
print(f'baseline: {nobs} inconsistent rows ({time.time()-t0:.0f}s)')
cand = collections.Counter()
for c in certs: cand.update(c)
cand = {a: n for a, n in cand.items() if a not in BAD}
print(f'candidate sacrifice atoms ({len(cand)}): '
      f'{sorted(((a, len(L.atom2eq.get(a,{}))) for a in cand), key=lambda t: t[1])[:24]}')

def cost(drop):
    lost = set()
    for a in drop: lost |= set(L.atom2eq.get(a, {}))
    return len(lost - GAIN), len(GAIN) - len(lost - GAIN)

# greedy: repeatedly drop the atom that removes the most obstructions per equation
drop = set()
for step in range(8):
    nobs, certs = analyse(drop)
    c, net = cost(drop)
    print(f'\nstep {step}: drop={sorted(drop)} cost={c} eqs, net={net:+d}, obstructions={nobs}')
    if nobs == 0:
        print(f'*** CONSISTENT with sacrifice {sorted(drop)}: net gain {net:+d} equations')
        break
    pool = collections.Counter()
    for cc in certs: pool.update(cc)
    opts = [a for a in pool if a not in BAD and a not in drop]
    opts.sort(key=lambda a: (len(set(L.atom2eq.get(a, {})) - GAIN), -pool[a]))
    best = None
    for a in opts[:10]:
        n2, _ = analyse(drop | {a}, want_cert=False)
        c2, net2 = cost(drop | {a})
        print(f'    try a{a:<6} (+{len(set(L.atom2eq.get(a,{}))-GAIN)} eqs) -> obstructions {n2}, net {net2:+d}')
        if best is None or (n2, c2) < best[0]: best = ((n2, c2), a)
    if best is None: break
    drop.add(best[1])
