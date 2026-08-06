"""CL: extract the minimal-support obstruction functionals and their escape routes."""
import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
vm = [x % P for x in v0]
av0 = L.all_atom_values(v0)
BAD = [a for a in range(L.NA) if a not in atom_out and av0[a]]

def grad_combo(weights):
    lam = collections.defaultdict(int)
    for a, wa in weights.items():
        if wa % P == 0: continue
        for w in L.avars[a]:
            lam[w] = (lam[w] + wa*ad.dpart(a, w, vm)) % P
    for t in reversed(ad.ORDER):
        lt = lam.get(t, 0) % P
        if not lt: continue
        a = definer[t]
        d = ad.dpart(a, t, vm)
        if d % P == 0: continue
        f = -lt*pow(d, -1, P) % P
        for w in L.avars[a]:
            if w == t: continue
            dw = ad.dpart(a, w, vm)
            if dw: lam[w] = (lam[w] + f*dw) % P
        lam[t] = 0
    return {u: lam[u] % P for u in ad.FREESET if lam.get(u, 0) % P}

D = json.load(open(os.path.join(HERE,'cl_closure_cols.json')))
cols = {int(u): {int(c): int(d) for c, d in m.items()} for u, m in D['cols'].items()}
U = sorted(cols)
rows = sorted(set().union(*[set(c) for c in cols.values()]) | set(BAD))
n, m = len(rows), len(U)
ri = {a: i for i, a in enumerate(rows)}
print(f'closure {n} x {m}')
A = [[0]*(m+1) for _ in range(n)]
for j, u in enumerate(U):
    for a, d in cols[u].items(): A[ri[a]][j] = d % P
for a in rows: A[ri[a]][m] = (-av0[a]) % P
Tr = [{i: 1} for i in range(n)]
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
        Ai = A[i]; A[i] = [(Ai[k]-f*Ar[k]) % P for k in range(m+1)]
        ti = Tr[i]
        for k, x in Trr.items():
            nv = (ti.get(k, 0)-f*x) % P
            if nv: ti[k] = nv
            elif k in ti: del ti[k]
    r += 1
print(f'rank {r}')
obs = [i for i in range(r, n) if A[i][m] % P and not any(A[i][j] for j in range(m))]
cands = sorted(obs, key=lambda i: len(Tr[i]))
out = []
for i in cands:
    y = {rows[k]: w for k, w in Tr[i].items() if w % P}
    g = grad_combo(y)
    val = sum(w*(av0[a] % P) for a, w in y.items()) % P
    nzsup = [a for a in y if av0[a]]
    print(f'\ncertificate |support|={len(y)}  y.r={val}')
    print(f'   escape free inputs: {sorted(g)}  (values ' +
          ', '.join(f'x_{u}={v0[u]}' for u in sorted(g)) + ')')
    print(f'   atoms in support with nonzero residue: {nzsup}')
    print(f'   support atoms: {sorted(y)}')
    out.append({'y': {str(a): str(w) for a, w in y.items()}, 'yr': str(val),
                'escapes': sorted(g)})
    if len(out) >= 3: break
json.dump(out, open(os.path.join(HERE,'cl_cert.json'),'w'))
print('\nsaved cl_cert.json')
