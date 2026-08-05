"""Step 18: close BOTH congruences on branch a -> residual {31670, 31672}. Measure."""
import sys
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

CODES, _ = H.load_equations()
H688 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
H1618 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002


def score(v, tag):
    nz = nz_checks(v); ng = nz_gates(v)
    live = [a for a in range(len(polys)) if evalpoly(polys[a], v) != 0]
    ff = H.evaluate(CODES, v, eqs_of(live))
    print(f'[{tag:36s}] nz={len(nz):3d} ng={len(ng)} FAIL={len(ff):4d} {sorted(set(nz+ng))}')
    return nz, ng, ff


v = H.load_assignment('quad/stateB1.json')   # 688 closed, {1618, 31672}
score(v, 'stateB1')

# --- close 1618: x_22649 and x_22152 shift together (keeps 2423) ---
k = (H1618 - v[30454]) % P
ripple(v, {22649: v[22649] + k, 22152: v[22152] + k})
print('x_30454 == H1618 (mod p):', (v[30454] - H1618) % P == 0, ' 2423 resid:', evalpoly(polys[2423], v))
score(v, 'x_22649,x_22152 += k')
d = v[24468] - H1618
if d % P == 0:
    ripple(v, {11436: d // P}); score(v, ' + x_11436 handle')
d2 = 12604395 * (v[22649] - v[29524])
if d2 % P == 0:
    ripple(v, {14768: d2 // P}); score(v, ' + x_14768 handle')

nz, ng, ff = score(v, 'BOTH CLOSED')
d = pickle.load(open('atoms.pkl', 'rb')); eq_terms = d['eq_terms']
live = sorted(a for a in range(len(polys)) if evalpoly(polys[a], v) != 0)
import collections
cnt = collections.Counter()
for i in ff:
    m, sq, tl = eq_terms[i]
    cnt[frozenset(a for c, a in tl if a in live)] += 1
for kk, n in cnt.most_common():
    print(f'   {n:3d} eqs involve {sorted(kk)}')
for a in live:
    print(f'   atom {a} in {len(atom2eq.get(a,[]))} eqs : {src[a][:170]}')
H.save_assignment(v, 'quad/stateC1.json')
print('saved quad/stateC1.json')
