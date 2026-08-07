"""WR step 18: realise the k=3 knob deformation with an INTEGRAL t by scaling the
uniform part.  d = c*1 + t with t supported on 3 wire coordinates, zeroing 3 of
the twelve a37694 rows.  Then repair and measure."""
import os, sys, itertools, math, json
from fractions import Fraction
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
rowsum = {e: sum(rows[e].values()) for e in RE}
BAD12 = [e for e in RE if rowsum[e]]
cand = sorted(set().union(*[set(rows[e]) for e in BAD12]))

DET = dict(W.F3)
for u in WIRE:
    if u in L.definer:
        DET[u] = L.definer[u]
FW = W.Frame(DET)
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
b2 = list(base); FW.fwd(b2)


def solveQ(coords, tgt, c):
    k = len(coords)
    M = [[Fraction(rows[e].get(j, 0)) for j in coords] + [Fraction(-c * rowsum[e])]
         for e in tgt]
    piv, r = [], 0
    for col in range(k):
        s = next((i for i in range(r, k) if M[i][col]), None)
        if s is None:
            return None
        M[r], M[s] = M[s], M[r]
        inv = M[r][col]
        M[r] = [x / inv for x in M[r]]
        for i in range(k):
            if i != r and M[i][col]:
                f = M[i][col]
                M[i] = [M[i][j] - f * M[r][j] for j in range(k + 1)]
        piv.append(col); r += 1
    return [M[i][k] for i in range(k)]


def broken(dv):
    return [e for e in RE if sum(cc * dv[j] for j, cc in rows[e].items() if dv[j])]


best = None
for k in (2, 3, 4):
    for coords in itertools.combinations(cand, k):
        for tgt in itertools.combinations(BAD12, k):
            t = solveQ(list(coords), list(tgt), 1)
            if t is None:
                continue
            den = 1
            for x in t:
                den = den * x.denominator // math.gcd(den, x.denominator)
            # c must be a multiple of den for t to be integral (t scales with c)
            d = [den] * N
            for i, j in enumerate(coords):
                d[j] += int(t[i] * den)
            bad = broken(d)
            if best is None or (len(bad), den) < (len(best[0]), best[3]):
                best = (bad, d, (coords, tgt), den)
bad, d, info, den = best
print(f'best integral knob deformation: {len(bad)} identity rows broken')
print(f'   coords {[WIRE[j] for j in info[0]]}, targets {info[1]}, scale c={den}')
print(f'   broken rows: {bad}')
print(f'   wire base value w = p + {den}; max |d| = {max(abs(x) for x in d)}')

v = list(b2)
for j, u in enumerate(WIRE):
    v[u] = P + d[j]
FW.fwd(v, rounds=10)
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
fail = L.failing_eqs(av)
print(f'MEASURED: score={L.NEQ-len(fail)} failing={len(fail)} nonzero={len(nz)} '
      f'{sorted(nz)[:25]}')
print(f'   failing: {sorted(fail)}')
T.save(v, os.path.join(HERE, 'wr_knob4.json'))
