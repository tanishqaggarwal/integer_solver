"""IP #10 -- the TRUE invariant: the least d such that  M x = d * rhs  has an integer solution.

The RREF particular solution's denominator is not necessarily minimal (the system has a
kernel).  The honest invariant is the smallest positive d with d*rhs in the integer column
span of M.  Test the divisors of the observed denominator D = 2458959 * p.
"""
import sys, os, json, time, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from zsolve import solve_int
from ip8 import build
from ip7 import load_raw
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
P = L.P

LAB = os.path.join(HERE, '..')
src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(LAB, 'best', 'new_instance_partial_39026.json')
v = load_raw(src)
print(f"=== {os.path.basename(src)}")
v, FAIL, used, M, rhs, nf = build(v)

cands = [1, 3, 819653, 2458959, P, 3 * P, 819653 * P, 2458959 * P]
print("least d with  M x = d*rhs  integer-solvable:")
for d in cands:
    t0 = time.time()
    x = solve_int(M, [d * r for r in rhs])
    lab = ('1' if d == 1 else ('P' if d == P else
           ('%d' % d if d < 10**9 else '%d*P' % (d // P))))
    print(f"   d = {lab:12s} : {'SOLVABLE' if x else 'no':9s} ({time.time()-t0:.0f}s)", flush=True)
    if x:
        print(f"   => the whole obstruction is a single divisibility by {lab}")
        break
