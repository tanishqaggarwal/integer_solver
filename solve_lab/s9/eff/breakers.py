"""Scan ALL 38,748 variables for moves that break one of the two surviving mod-P
congruences of the S13 lattice, and price each in extra sacrificed equations.

congruence 1 :  y(22230) + y(22231)                   (mod P)
congruence 2 :  y(22229) + 7376877 * y(35762)         (mod 7376877 * P)
"""
import pickle, sys, time, collections
import lib as L, model as MD

v0 = L.load(L.BEST24)
MD.BASEP = [MD.prim_val(a, v0) for a in range(L.NA)]
S13 = frozenset([2554, 6816, 8124, 8680, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125])
A = set(MD.confined_atoms(S13))
P = L.P
M2 = 7376877 * P
DELTA = 1

rows = []
t0 = time.time()
for x in range(L.NVARS):
    v, tou = MD.move(v0, {x: v0[x] + DELTA}, A)
    if not tou:
        continue
    d = {a: tou[a] - MD.BASEP[a] for a in tou}
    c1 = (d.get(22230, 0) + d.get(22231, 0)) % P
    c2 = (d.get(22229, 0) + 7376877 * d.get(35762, 0)) % M2
    if c1 == 0 and c2 == 0:
        continue
    extra_atoms = sorted(set(tou) - A)
    cost = sorted(L.eqs_of_atoms(extra_atoms) - S13)
    rows.append((len(cost), x, c1 != 0, c2 != 0, extra_atoms, cost))
    if x % 8000 == 0:
        print(f'  {x}/{L.NVARS} {time.time()-t0:.0f}s', file=sys.stderr)
rows.sort()
print(f'variables that move a congruence: {len(rows)}   ({time.time()-t0:.0f}s)')
for cost, x, b1, b2, ea, ce in rows[:60]:
    print(f'  x_{x:<6d} cost={cost:<3d} breaks C1={b1} C2={b2}  extra_atoms={ea}  eqs={ce[:20]}')
pickle.dump(rows, open('breakers.pkl', 'wb'))
