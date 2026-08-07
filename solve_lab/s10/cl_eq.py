"""CL: the EQUATION-level mod-p system.  Equations, not atoms, must vanish, so atoms are
allowed to cancel inside an equation.  This is strictly weaker than the atom closure."""
import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0 = L.all_atom_values(v0)

D = json.load(open(os.path.join(HERE,'cl_closure_cols.json')))
cols = {int(u): {int(c): int(d) for c, d in m.items()} for u, m in D['cols'].items()}
U = sorted(cols); m = len(U)
atoms = sorted(set().union(*[set(c) for c in cols.values()]))
aset = set(atoms)
EQ = sorted(set().union(*[set(L.atom2eq.get(a, {})) for a in atoms]))
print(f'{len(atoms)} moving atoms -> {len(EQ)} equations x {m} free inputs')

n = len(EQ)
A = [[0]*(m+1) for _ in range(n)]
for i, e in enumerate(EQ):
    mult, sq, co = L.eq_atoms[e]
    s = 0
    for a, c in co.items():
        s += c*av0[a]
        if a in aset:
            for j, u in enumerate(U):
                d = cols[u].get(a)
                if d: A[i][j] = (A[i][j] + c*d) % P
    A[i][m] = (-s) % P
nzr = sum(1 for i in range(n) if any(A[i][:m]) or A[i][m])
print(f'rows with any content: {nzr}')

Tr = [{i: 1} for i in range(n)]
t0 = time.time()
r = 0
for c in range(m):
    pr = next((i for i in range(r, n) if A[i][c]), None)
    if pr is None: continue
    A[r], A[pr] = A[pr], A[r]; Tr[r], Tr[pr] = Tr[pr], Tr[r]
    inv = pow(A[r][c], -1, P)
    if inv != 1:
        A[r] = [x*inv % P for x in A[r]]; Tr[r] = {k: x*inv % P for k, x in Tr[r].items()}
    Ar, Trr = A[r], Tr[r]
    for i in range(r+1, n):
        f = A[i][c]
        if not f: continue
        A[i] = [(A[i][k]-f*Ar[k]) % P for k in range(m+1)]
        ti = Tr[i]
        for k, x in Trr.items():
            nv = (ti.get(k, 0)-f*x) % P
            if nv: ti[k] = nv
            elif k in ti: del ti[k]
    r += 1
    if r % 100 == 0: print(f'  pivot {r} ({time.time()-t0:.0f}s)', flush=True)
print(f'rank {r} of {m} cols, {n} rows  ({time.time()-t0:.0f}s)')
obs = [i for i in range(r, n) if A[i][m] % P and not any(A[i][j] for j in range(m))]
print(f'INCONSISTENT equation-combinations: {len(obs)}')
if not obs:
    print('*** EQUATION-LEVEL SYSTEM IS CONSISTENT -- full repair exists mod p!')
    sol = [0]*m
    for i in range(len(A)-1, -1, -1): pass
    # back substitute using the echelon rows
    piv = []
    for i in range(r):
        c = next(j for j in range(m) if A[i][j])
        piv.append((i, c))
    for i, c in reversed(piv):
        s = A[i][m]
        for j in range(c+1, m):
            if A[i][j] and sol[j]: s = (s - A[i][j]*sol[j]) % P
        sol[c] = s % P
    json.dump({str(U[j]): str(sol[j]) for j in range(m) if sol[j]},
              open(os.path.join(HERE, 'cl_eqdelta.json'), 'w'))
    print(f'saved cl_eqdelta.json  ({sum(1 for x in sol if x)} inputs move)')
else:
    cnt = collections.Counter()
    minsup = None
    for i in obs:
        y = {EQ[k]: w for k, w in Tr[i].items() if w % P}
        cnt.update(y)
        if minsup is None or len(y) < len(minsup): minsup = y
    print(f'  smallest obstruction combines {len(minsup)} equations: {sorted(minsup)}')
    print(f'  equations appearing in obstructions: {len(cnt)}')
    print(f'  most frequent: {cnt.most_common(20)}')
    json.dump({'obs': len(obs), 'min_support': [str(k) for k in sorted(minsup)],
               'freq': {str(k): c for k, c in cnt.most_common(60)}},
              open(os.path.join(HERE, 'cl_eqobs.json'), 'w'))
