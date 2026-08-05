"""Step 21: 2-D bit lattice search over (m30, m19) and a cost test for activating bits."""
import sys, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

CODES, _ = H.load_equations()
dd = pickle.load(open('atoms.pkl', 'rb')); eq_terms = dd['eq_terms']
D19 = 77950674053801246186103680324559083425863605288488642042743492031016466801949
D30 = 33784113124339387961669114809489170549755964517992654161423235401600078441306
t19 = 77995533598565078661335603578384652517982926814434273088397111539457287830836
t30 = 69985052393467842547745661261097129016803957409171160345035789691482532567532
res = pickle.load(open('quad/sens.pkl', 'rb'))
both = [t for t, r in res.items() if r[0] == D19 and r[1] == D30]
only19 = [t for t, r in res.items() if r[0] == D19 and r[1] == 0]
only30 = [t for t, r in res.items() if r[0] == 0 and r[1] == D30]
print(f'bits: both={len(both)} only19={len(only19)} only30={len(only30)}')
NA, NC, NB = len(both), len(only19), len(only30)
M19MAX, M30MAX = NA + NC, NA + NB
print(f'reachable m19 in [0,{M19MAX}], m30 in [0,{M30MAX}] (coupled through the {NA} shared bits)')

EQS = [(4833, 34, -11, 1), (4944, 1, 6, 0), (5348, -23, 4, 0), (9344, -16, -2, 0),
       (10406, 1, 0, 0), (11574, 29, -20, 1), (12321, 23, 11, 0), (19708, 0, 7, 0),
       (20927, 0, 1, 0), (21972, 21, 7, 0), (27514, 40, 31, 0), (38014, 35, 34, 0)]

print('\n--- can any failing equation be zeroed by a bit choice? ---')
hits = []
for i, c1, c2, sq in EQS:
    found = []
    for m30 in range(0, M30MAX + 1):
        A = (t30 - m30 * D30) % P
        for m19 in range(0, M19MAX + 1):
            B = (t19 - m19 * D19) % P
            if (c1 * A + c2 * B) % P == 0:
                # feasibility of (m19,m30) with the coupling
                na_lo = max(0, m19 - NC, m30 - NB)
                na_hi = min(NA, m19, m30)
                if na_lo <= na_hi:
                    found.append((m30, m19))
    print(f'  eq {i:6d} (c1={c1},c2={c2}): {len(found)} bit-configs zero it {found[:4]}')
    if found: hits.append((i, found))

print('\n--- can a whole defect vanish? ---')
for nm, t, D, mx in (('A (31670)', t30, D30, M30MAX), ('B (31672)', t19, D19, M19MAX)):
    sol = [m for m in range(0, mx + 1) if (t - m * D) % P == 0]
    print(f'  {nm} == 0 mod p for m in [0,{mx}]: {sol}')
    # also with an extra p-multiple of the pinned var (no effect mod p)
print('\n--- cost of activating one bit at stateC1 ---')
v0 = H.load_assignment('quad/stateC1.json')
base_nz = set(nz_checks(v0))
for b in (both[:3] + only19[:2] + only30[:2]):
    v = list(v0); ripple(v, {b: 1 - v0[b]})
    nz = set(nz_checks(v))
    live = [a for a in range(len(polys)) if evalpoly(polys[a], v) != 0]
    ff = H.evaluate(CODES, v, eqs_of(live))
    print(f'  bit x_{b}: nz={len(nz)} FAIL={len(ff)} new={sorted(nz-base_nz)[:10]}')
