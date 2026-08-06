"""WR step 19: search for a wire deformation d = (1-p)*1 + t (wire base value 1,
so the handles stay unquantised) with t INTEGRAL, supported on three wire
coordinates, minimising the broken rows among the 219 identity rows PLUS the
a39417 row (eq 11915)."""
import os, sys, itertools, math, json
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_rows as R
P = ad.P
WIRE, widx, rows, RE = R.WIRE, R.widx, R.rows, R.RE
N = len(WIRE)

# augment with the a39417 row (its equation 11915 vanishes iff the form does)
f39417 = {}
for m, c in L.polys[39417].items():
    if len(m) == 1 and m[0] in widx:
        f39417[widx[m[0]]] = f39417.get(widx[m[0]], 0) + c
ROWS = dict(rows)
ROWS[-1] = f39417                        # pseudo-row for eq 11915
ALL = sorted(ROWS)
rowsum = {e: sum(ROWS[e].values()) for e in ALL}
BAD = [e for e in ALL if rowsum[e]]
print(f'rows broken by the uniform shift: {len(BAD)} -> {BAD}')
cand = sorted(set().union(*[set(ROWS[e]) for e in BAD]))
print(f'candidate coordinates: {len(cand)}')
C = 1 - P


def det3(a):
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


def broken(dv):
    return [e for e in ALL
            if sum(cc * dv[j] for j, cc in ROWS[e].items() if dv[j])]


best = None
tested = solvable = 0
for coords in itertools.combinations(cand, 3):
    cols = list(coords)
    for tgt in itertools.combinations(BAD, 3):
        A = [[ROWS[e].get(j, 0) for j in cols] for e in tgt]
        D = det3(A)
        if D == 0:
            continue
        tested += 1
        b = [-rowsum[e] for e in tgt]          # per unit c
        ts = []
        ok = True
        for i in range(3):
            Ai = [row[:] for row in A]
            for r in range(3):
                Ai[r][i] = b[r]
            Di = det3(Ai)
            num = C * Di
            if num % D:
                ok = False
                break
            ts.append(num // D)
        if not ok:
            continue
        solvable += 1
        d = [C] * N
        for i, j in enumerate(cols):
            d[j] += ts[i]
        bad = broken(d)
        if best is None or len(bad) < len(best[0]):
            best = (bad, d, cols, tgt)
            print(f'  new best {len(bad)} broken  coords '
                  f'{[WIRE[j] for j in cols]} targets {tgt}', flush=True)
print(f'\ntested {tested} systems, {solvable} with an integral solution at wire base 1')
if best:
    bad, d, cols, tgt = best
    print(f'BEST: {len(bad)} rows broken -> {bad}')
    print(f'   knob coordinates {[WIRE[j] for j in cols]}, values '
          f'{[len(str(abs(d[j]))) for j in cols]} digits')
    json.dump({'d': [str(x) for x in d], 'bad': bad,
               'coords': [WIRE[j] for j in cols]},
              open(os.path.join(HERE, 'wr_knob5.json'), 'w'))
