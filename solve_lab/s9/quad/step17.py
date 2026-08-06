"""Step 17: full reconstruction with the quadrant bits FROZEN."""
import sys
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
import driver
from driver import *
driver.CODES, _ = H.load_equations()
driver.FREEZE = {2081, 24601}

lbl = sys.argv[1]
seeds = {'a': {2081: 0}, 'b': {24601: 0}, 'ab': {2081: 0, 24601: 0}}[lbl]
v = H.load_assignment(BEST)
ripple(v, seeds)
driver.score(v, f'branch {lbl} raw')
v, hist = driver.repair(v, rounds=40)
nz, ng, ff = driver.score(v, f'branch {lbl} repaired')
print('history:', hist)
for a in sorted(set(nz + ng)):
    print(f'   atom {a}: {src[a][:180]}')
H.save_assignment(v, f'quad/frz_{lbl}.json')
print('saved quad/frz_%s.json' % lbl)
