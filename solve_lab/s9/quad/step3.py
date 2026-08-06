"""Step 3: repair the 'needs ZERO' free inputs on the x_2081=0 branch, then iterate."""
import sys, time, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

CODES, _ = H.load_equations()
v0 = H.load_assignment(BEST)
v = list(v0)
ripple(v, {2081: 0})


def score(v, tag):
    nz = nz_checks(v); ng = nz_gates(v)
    live = [a for a in range(len(polys)) if evalpoly(polys[a], v) != 0]
    ff = H.evaluate(CODES, v, eqs_of(live))
    print(f'[{tag}] nzcheck={len(nz)} nzgate={len(ng)} liveraw={len(live)} FAIL={len(ff)}')
    return nz, ng, ff


score(v, 'branch a')
for u in (24548, 14623, 14853, 31339):
    ripple(v, {u: 0})
    score(v, f'after x_{u}:=0')

print('\nremaining residual atoms:')
nz, ng, ff = score(v, 'final')
for a in nz + ng:
    R = resid_poly.get(a, polys[a])
    print(f'  atom {a}: {src[a][:150]}')
H.save_assignment(v, 'quad/stateA1.json')
print('saved quad/stateA1.json')
