"""S12 step 4: closure diagnostic -- WHICH rows are inconsistent, and what is
the ceiling of the 39009 frame anyway?
"""
import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE); FORBID = {2081, 4287}
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
vm0 = [x % P for x in v0]
av0 = L.all_atom_values(v0)
BAD = [21617, 29539]
NZ = [a for a in range(L.NA) if av0[a]]
CHECKS = sorted(a for a in range(L.NA) if a not in atom_out)
print('nonzero atoms', NZ)
for a in NZ:
    e = sorted(L.atom2eq.get(a, {}))
    print(f'  a{a}: {len(e)} eqs {e}')
allf = set(L.failing_eqs(av0))
for sub in ([21617,29539], [37662,40826], [21617], [29539], [37662], [40826]):
    rest = [a for a in NZ if a not in sub]
    e = sorted(L.eqs_of_atoms(rest))
    print(f'  if {sub} were zeroed: remaining eqs touched {len(e)} -> ceiling score <= {L.NEQ-0 if not e else L.NEQ-len(e)}')

# ---- closure ---------------------------------------------------------------
def closure_of(v, extra_bad=()):
    vm = [x % P for x in v]
    gc = {}
    def gr(c):
        if c not in gc: gc[c] = set(ad.grad(c, vm)) - FORBID
        return gc[c]
    U, rows, cols = set(), set(BAD) | set(extra_bad), {}
    for it in range(9):
        newU = set()
        for c in rows: newU |= gr(c)
        newU -= U; U |= newU
        for u in sorted(newU): cols[u] = jac_column(u, v, vm, CHECKS)
        nr = set(BAD) | set(extra_bad)
        for u in U: nr |= set(cols[u])
        grew = nr - rows; rows = nr
        if not newU and not grew: break
    return sorted(rows), sorted(U), cols

t0 = time.time()
rows, Us, cols = closure_of(v0)
print(f'\nclosure {len(rows)} rows x {len(Us)} cols  ({time.time()-t0:.0f}s)')
print(f'  nonzero-residual rows inside the closure: {[r for r in rows if av0[r]]}')

ri = {c:i for i,c in enumerate(rows)}
n, m = len(rows), len(Us)
M = [[0]*(m+1) for _ in rows]
for j,u in enumerate(Us):
    for c,d in cols[u].items():
        if c in ri: M[ri[c]][j] = d % P
for c in rows: M[ri[c]][m] = (-av0[c]) % P
# track row provenance through elimination
prov = [{i:1} for i in range(n)]
r_ = 0; piv = []
for j in range(m):
    k = next((i for i in range(r_, n) if M[i][j]), None)
    if k is None: continue
    M[r_],M[k] = M[k],M[r_]; prov[r_],prov[k] = prov[k],prov[r_]
    inv = pow(M[r_][j], -1, P)
    M[r_] = [x*inv % P for x in M[r_]]
    prov[r_] = {a:b*inv % P for a,b in prov[r_].items()}
    for i in range(n):
        if i != r_ and M[i][j]:
            f = M[i][j]
            M[i] = [(a2-f*b2) % P for a2,b2 in zip(M[i], M[r_])]
            pp = dict(prov[i])
            for a,b in prov[r_].items(): pp[a] = (pp.get(a,0)-f*b) % P
            prov[i] = {a:b for a,b in pp.items() if b}
    piv.append(j); r_ += 1
badrows = [i for i in range(r_, n) if M[i][m]]
print(f'  rank {r_} / cols {m}  kernel {m-r_}  INCONSISTENT rows {len(badrows)}')
print('\n  the inconsistency functionals (checks involved, weight-sorted):')
for i in badrows[:12]:
    sup = sorted(prov[i], key=lambda k: -abs(prov[i][k]))
    ch = [rows[k] for k in sup]
    print(f'    #{i}: support {len(ch)} checks; nonzero-residual members '
          f'{[rows[k] for k in sup if av0[rows[k]]]}')
# common structure
allsup = collections.Counter()
for i in badrows:
    for k in prov[i]: allsup[rows[k]] += 1
print(f'  checks appearing in >= half the inconsistency functionals: '
      f'{[c for c,n2 in allsup.items() if n2 >= len(badrows)/2][:40]}')
json.dump({'rows':rows,'Us':Us,'badrows':len(badrows),'rank':r_}, open(os.path.join(HERE,'ac_diag.json'),'w'))
