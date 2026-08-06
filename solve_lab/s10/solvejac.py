"""S10 step 102: solve J.delta == -r  (mod p) over the 572x134 closure."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
J = json.load(open(os.path.join(HERE, 'jac.json')))
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
U = J['U']; rows = J['rows']
cols = {int(u): {int(c): int(d) for c, d in m.items()} for u, m in J['cols'].items()}
resid = {int(c): int(d) for c, d in J['resid'].items()}
ri = {c: i for i, c in enumerate(rows)}
M = [[0] * (len(U) + 1) for _ in rows]
for j, u in enumerate(U):
    for c, d in cols[u].items():
        M[ri[c]][j] = d % P
for c in rows:
    M[ri[c]][len(U)] = (-resid.get(c, 0)) % P
n, m = len(rows), len(U)
print(f'system {n} x {m} over F_p; rhs nonzero rows: '
      f'{[c for c in rows if resid.get(c, 0)]}')

piv, r_ = [], 0
for j in range(m):
    k = next((i for i in range(r_, n) if M[i][j]), None)
    if k is None: continue
    M[r_], M[k] = M[k], M[r_]
    inv = pow(M[r_][j], -1, P)
    M[r_] = [x * inv % P for x in M[r_]]
    for i in range(n):
        if i != r_ and M[i][j]:
            f = M[i][j]
            M[i] = [(a - f * b) % P for a, b in zip(M[i], M[r_])]
    piv.append(j); r_ += 1
rank = r_
incons = [i for i in range(rank, n) if M[i][m]]
print(f'rank(J) = {rank};  inconsistent rows after elimination: {len(incons)}')
if incons:
    print('  -> the linearised repair is INCONSISTENT mod p')
    for i in incons[:5]:
        print(f'     witness row {i}')
else:
    print('  -> CONSISTENT.  extracting a solution')
    delta = [0] * m
    for i, j in enumerate(piv):
        delta[j] = M[i][m]
    nzu = [(U[j], delta[j]) for j in range(m) if delta[j]]
    print(f'  solution moves {len(nzu)} free inputs')
    json.dump({str(u): str(d) for u, d in nzu},
              open(os.path.join(HERE, 'delta.json'), 'w'))
    print('  saved delta.json')
    print(f'  kernel dimension {m - rank}')
