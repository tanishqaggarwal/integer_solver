"""CL: price the obstruction certificate, and test each congruence separately."""
import os, sys, json, collections
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

CERT = [1436, 3576, 3578, 7930, 7932, 15456, 15462, 21617, 21619, 40065, 41507]
print('=== the 11-atom obstruction certificate ===')
E21 = set(L.atom2eq.get(21617, {}))
for a in CERT:
    eqs = set(L.atom2eq.get(a, {}))
    print(f'  a{a:<6} eqs={len(eqs):<3} overlap_with_a21617={len(eqs & E21):<3}  {L.atom_src[a][:120]}')
print(f'\na21617 eqs = {sorted(E21)} ({len(E21)}) ; shadow a37662 eqs {sorted(L.atom2eq.get(37662,{}))}')
gain = len(E21 | set(L.atom2eq.get(37662, {})))
print(f'zeroing a21617 (+shadow a37662) recovers {gain} equations')
print('\nprice of breaking each certificate member instead:')
best = None
for a in CERT:
    if a == 21617: continue
    eqs = set(L.atom2eq.get(a, {}))
    net = gain - len(eqs - (E21 | set(L.atom2eq.get(37662, {}))))
    print(f'  break a{a:<6}: costs {len(eqs)} eqs ({len(eqs - E21)} outside a21617) -> net {net:+d}')
    if best is None or net > best[0]: best = (net, a)
print(f'  BEST net: {best[0]:+d} by breaking a{best[1]}')

# ---------- consistency of each congruence alone ----------
D = json.load(open(os.path.join(HERE,'cl_closure_cols.json')))
cols = {int(u): {int(c): int(d) for c, d in m.items()} for u, m in D['cols'].items()}
U = sorted(cols)
allrows = sorted(set().union(*[set(c) for c in cols.values()]) | set(BAD))

def solve(drop):
    rows = [a for a in allrows if a not in drop]
    n, m = len(rows), len(U)
    ri = {a: i for i, a in enumerate(rows)}
    A = [[0]*(m+1) for _ in range(n)]
    for j, u in enumerate(U):
        for a, d in cols[u].items():
            if a in ri: A[ri[a]][j] = d % P
    for a in rows: A[ri[a]][m] = (-av0[a]) % P
    r = 0; piv = []
    for c in range(m):
        pr = next((i for i in range(r, n) if A[i][c]), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]
        inv = pow(A[r][c], -1, P)
        if inv != 1: A[r] = [x*inv % P for x in A[r]]
        Ar = A[r]
        for i in range(r+1, n):
            f = A[i][c]
            if f: A[i] = [(A[i][k]-f*Ar[k]) % P for k in range(m+1)]
        piv.append(c); r += 1
    bad = sum(1 for i in range(r, n) if A[i][m] % P and not any(A[i][j] for j in range(m)))
    return r, n, bad, A, piv, rows

for drop, tag in [(set(), 'FULL: both congruences'),
                  ({21617, 37662, 21619}, 'a29539 only (a21617 released)'),
                  ({29539, 40826}, 'a21617 only (a29539 released)'),
                  ({21617, 37662, 21619, 3576}, 'a29539 + a21617 released except a3576 broken')]:
    r, n, bad, A, piv, rows = solve(drop)
    print(f'\n[{tag}] rows {n} rank {r}  inconsistent {bad}  -> '
          f'{"CONSISTENT" if bad == 0 else "INCONSISTENT"}')
    if bad == 0 and drop:
        m = len(U)
        sol = [0]*m
        for i in range(len(piv)-1, -1, -1):
            c = piv[i]; s = A[i][m]
            for j in range(c+1, m):
                if A[i][j] and sol[j]: s = (s - A[i][j]*sol[j]) % P
            sol[c] = s % P
        nz = {U[j]: sol[j] for j in range(m) if sol[j]}
        print(f'   solution moves {len(nz)} free inputs')
        json.dump({str(k): str(x) for k, x in nz.items()},
                  open(os.path.join(HERE, f'cl_delta_{"_".join(map(str,sorted(drop)))}.json'), 'w'))
