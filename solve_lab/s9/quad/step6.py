"""Step 6: dependency-ordered repair on the x_2081=0 branch."""
import sys
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

CODES, _ = H.load_equations()
H688 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
H1618 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002

v = H.load_assignment('quad/stateA1.json')


def score(v, tag, show=False):
    nz = nz_checks(v); ng = nz_gates(v)
    live = [a for a in range(len(polys)) if evalpoly(polys[a], v) != 0]
    ff = H.evaluate(CODES, v, eqs_of(live))
    print(f'[{tag:28s}] nzcheck={len(nz):3d} nzgate={len(ng)} FAIL={len(ff):4d}  {(nz+ng)[:20]}')
    return nz, ng, ff


score(v, 'start')
for u, val in [(6418, 0), (12553, 0), (22162, 0), (30213, 0)]:
    ripple(v, {u: val})
    score(v, f'x_{u}:=0')

print('\n--- now the 688 / 1618 pins ---')
print('x_18956 =', v[18956] == v[16742] + P * v[22820], '(= x_16742 + p*x_22820)')
print('x_24468 =', v[24468] == v[12186] + P * 12354891 * v[14393], '(= x_12186 + p*12354891*x_14393)')
print('x_16742 mod p == H688 mod p ?', (v[16742] - H688) % P == 0)
print('x_12186 mod p == H1618 mod p ?', (v[12186] - H1618) % P == 0)
print('x_19083 =', v[19083], 'x_16742 =', v[16742], 'equal:', v[19083] == v[16742])
print('x_6361 def', definer.get(6361), 'x_23758 def', definer.get(23758))
print('x_33787 def', definer.get(33787), ' x_38100 val==p:', v[38100] == P)
print('atom 26731 resid:', evalpoly(polys[26731], v))

H.save_assignment(v, 'quad/stateA2.json')
print('saved quad/stateA2.json')
