"""Step 11: close atom 1618 via x_22649 (mod-p solve) + its pin partner x_22152."""
import sys
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

CODES, _ = H.load_equations()
H688 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
H1618 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
v0 = H.load_assignment('quad/stateA2.json')


def score(v, tag):
    nz = nz_checks(v); ng = nz_gates(v)
    live = [a for a in range(len(polys)) if evalpoly(polys[a], v) != 0]
    ff = H.evaluate(CODES, v, eqs_of(live))
    print(f'[{tag:38s}] nzcheck={len(nz):3d} nzgate={len(ng)} FAIL={len(ff):4d}  {sorted(set(nz+ng))[:22]}')
    return nz, ng, ff


score(v0, 'stateA2')

# exact coefficient of x_22649 on x_30454
v = list(v0); ripple(v, {22649: v0[22649] + 1})
c30 = v[30454] - v0[30454]
print('d x_30454 / d x_22649 =', c30, ' invertible mod p:', c30 % P != 0)
need = (H1618 - v0[30454]) % P
k = need * pow(c30 % P, -1, P) % P
print('k =', k)

v = list(v0)
ripple(v, {22649: v0[22649] + k})
print('x_30454 == H1618 mod p:', (v[30454] - H1618) % P == 0)
score(v, 'x_22649 += k')

# close the 1618 handles: x_24468 - H1618 = x_32989 = p*x_11436 ; x_12186 = x_30454 + p*x_33612
d = v[24468] - H1618
print('x_24468 - H1618 divisible by p:', d % P == 0)
if d % P == 0:
    ripple(v, {11436: d // P})
    score(v, 'x_11436 handle set')

# close the 2423 pin: need x_22649 == x_29524 (mod p);  x_29524 == x_22152 (free)
print('x_29524 =', v[29524] == v[22152], '(== x_22152)')
ripple(v, {22152: v[22649]})
score(v, 'x_22152 := x_22649')
d2 = 12604395 * (v[22649] - v[29524])
print('2423 numerator div by p:', d2 % P == 0, ' resid 2423 =', evalpoly(polys[2423], v))
if d2 % P == 0 and d2 != 0:
    ripple(v, {14768: d2 // P})
    score(v, 'x_14768 handle')

nz, ng, ff = score(v, 'FINAL')
for a in sorted(set(nz + ng)):
    print(f'   atom {a}: {src[a][:140]}')
H.save_assignment(v, 'quad/stateA3.json')
print('saved quad/stateA3.json')
