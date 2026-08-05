"""Step 15: automated reconstruction on each quadrant branch."""
import sys
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
import driver
from driver import *
driver.CODES, _ = H.load_equations()

lbl = sys.argv[1]
seeds = {'a': {2081: 0}, 'b': {24601: 0}, 'ab': {2081: 0, 24601: 0}, 'base': {}}[lbl]
v = H.load_assignment(BEST)
ripple(v, seeds)
driver.score(v, f'branch {lbl} raw')
v, hist = driver.repair(v, rounds=30)
nz, ng, ff = driver.score(v, f'branch {lbl} repaired')
print('history:', hist)
for a in sorted(set(nz + ng)):
    print(f'   atom {a}: {src[a][:170]}')
H.save_assignment(v, f'quad/auto_{lbl}.json')
print('saved quad/auto_%s.json' % lbl)
