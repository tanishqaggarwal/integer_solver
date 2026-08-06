"""WR step 20: realise + measure the best integral 3-knob deformation at wire
base 1 (11 identity-model rows broken, a39417 row saved)."""
import os, sys, itertools, math
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W
import wr_rows as R
P = ad.P
WIRE, widx, rows, RE = R.WIRE, R.widx, R.rows, R.RE
N = len(WIRE)
f39417 = {}
for m, c in L.polys[39417].items():
    if len(m) == 1 and m[0] in widx:
        f39417[widx[m[0]]] = f39417.get(widx[m[0]], 0) + c
ROWS = dict(rows); ROWS[-1] = f39417
ALL = sorted(ROWS)
rowsum = {e: sum(ROWS[e].values()) for e in ALL}
C = 1 - P
cols = [widx[u] for u in (12752, 13720, 18306)]
tgt = [-1, 12594, 25313]


def det3(a):
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


A = [[ROWS[e].get(j, 0) for j in cols] for e in tgt]
D = det3(A)
b = [-rowsum[e] for e in tgt]
ts = []
for i in range(3):
    Ai = [row[:] for row in A]
    for r in range(3):
        Ai[r][i] = b[r]
    ts.append(C * det3(Ai) // D)
d = [C] * N
for i, j in enumerate(cols):
    d[j] += ts[i]
bad = [e for e in ALL if sum(cc * d[j] for j, cc in ROWS[e].items() if d[j])]
print(f'identity-model rows broken: {len(bad)} -> {bad}')
print(f'wire values: base 1; knobs x_12752={len(str(abs(P+d[cols[0]])))} digits, '
      f'x_13720={len(str(abs(P+d[cols[1]])))}, x_18306={len(str(abs(P+d[cols[2]])))}')

DET = dict(W.F3)
for u in WIRE:
    if u in L.definer:
        DET[u] = L.definer[u]
FW = W.Frame(DET)
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
b2 = list(base); FW.fwd(b2)
v = list(b2)
for j, u in enumerate(WIRE):
    v[u] = P + d[j]
FW.fwd(v, rounds=10)
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
fail = L.failing_eqs(av)
print(f'MEASURED: score={L.NEQ-len(fail)} failing={len(fail)} nonzero={len(nz)}')
print(f'   nonzero atoms: {sorted(nz)}')
print(f'   failing eqs  : {sorted(fail)}')
T.save(v, os.path.join(HERE, 'wr_knob6.json'))
import wr_engine as E
E.Engine(FW, forbid=set(WIRE)).run(v, 'knob6', budget=1500, save=True)
