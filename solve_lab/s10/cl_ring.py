"""CL: the exact affine ring.  Restrict every cluster gadget + pin to the handful of
free inputs that carry the broadcast constants, and solve the little system over F_p."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
vm = [x % P for x in v0]
av0 = L.all_atom_values(v0)

K1 = v0[33462] % P
K2 = v0[22152] % P
NODES = [8778, 14623, 16742, 31339, 33462,      # K1 free class
         24548, 12553,                          # the a7930 chain
         14853, 22152, 22649,                   # K2 free class
         6418, 17325, 9413, 8731]               # other live free inputs in the cones
ATOMS = [(31672,'pin  x_33462 == K1'), (3578,'pin  x_12553 == C3'),
         (7930, 'gad  x_24548 == x_25442'), (21617,'gad  x_14623 == x_27522'),
         (33796,'gad  x_31339 == x_6858'), (26731,'gad  x_16742 == x_19083'),
         (33929,'gad  x_8778  == x_16144'),
         (31670,'pin  x_22152 == K2'), (29539,'gad  x_14853 == x_1308'),
         (2423, 'gad  x_22649 == x_29524'), (3576,'pin  x_6418  == C4')]

print(f'K1 mod p = {K1}\nK2 mod p = {K2}\n')
rowsdata = []
for a, desc in ATOMS:
    g = ad.grad(a, vm)
    full = len(g)
    row = {u: g.get(u, 0) % P for u in NODES if g.get(u, 0) % P}
    outside = {u: d for u, d in g.items() if u not in NODES and u not in (2081, 4287)}
    r = av0[a] % P
    rowsdata.append((a, desc, row, r, outside))
    print(f'a{a:<6} {desc:<28} resid={r!=0}  support {full} '
          f'(inside {len(row)}, outside {len(outside)})')
    print(f'      ' + '  '.join(f'x_{u}:{str(d)[:10]}..' for u, d in sorted(row.items())))
    if outside:
        print(f'      outside free inputs: {sorted(outside)[:16]}')

# ---- solve the little system over F_p:  for each atom, resid + sum d_u * delta_u = 0
cols = sorted(set().union(*[set(rd[2]) for rd in rowsdata]))
print(f'\nlittle system: {len(rowsdata)} rows x {len(cols)} cols  (cols {cols})')
M = [[rd[2].get(u, 0) for u in cols] + [(-rd[3]) % P] for rd in rowsdata]
lab = [rd[0] for rd in rowsdata]
n, m = len(M), len(cols)
Tr = [{i: 1} for i in range(n)]
r = 0
for c in range(m):
    pr = next((i for i in range(r, n) if M[i][c]), None)
    if pr is None: continue
    M[r], M[pr] = M[pr], M[r]; Tr[r], Tr[pr] = Tr[pr], Tr[r]
    inv = pow(M[r][c], -1, P)
    M[r] = [x*inv % P for x in M[r]]; Tr[r] = {k: x*inv % P for k, x in Tr[r].items()}
    for i in range(n):
        if i != r and M[i][c]:
            f = M[i][c]
            M[i] = [(M[i][k]-f*M[r][k]) % P for k in range(m+1)]
            for k, x in Tr[r].items():
                nv = (Tr[i].get(k, 0)-f*x) % P
                if nv: Tr[i][k] = nv
                elif k in Tr[i]: del Tr[i][k]
    r += 1
print(f'rank {r} of {m} columns / {n} rows')
for i in range(r, n):
    if M[i][m]:
        y = {lab[k]: w for k, w in Tr[i].items() if w}
        print(f'\n*** CONTRADICTION ROW: 0 = {M[i][m]}')
        print(f'    combination: ' + ' + '.join(f'{str(w)[:12]}..*a{a}' for a, w in y.items()))
        print(f'    i.e. these atoms cannot all vanish while the others hold')
    elif not any(M[i][:m]):
        y = {lab[k]: w for k, w in Tr[i].items() if w}
        print(f'\n    dependent row (consistent): ' + ' + '.join(f'a{a}' for a in y))
