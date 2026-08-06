"""Step 7: the two remaining congruences  x_16742 == H688 (p)  and  x_12186 == H1618 (p)."""
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
    print(f'[{tag:34s}] nzcheck={len(nz):3d} nzgate={len(ng)} FAIL={len(ff):4d}  {(nz+ng)[:22]}')
    return nz, ng, ff


score(v0, 'stateA2')

print('\n=== experiment: x_16742 := H688 ===')
v = list(v0)
ripple(v, {16742: H688})
print('  x_19083 now', v[19083] == H688, ' x_18956==H688:', v[18956] == H688)
nz, ng, ff = score(v, 'x_16742:=H688')

print('\n=== experiment: also fix 26731 handle x_33787 ===')
d = 6788513 * (v[16742] - v[19083])
print('  6788513*(x_16742-x_19083) =', d, ' divisible by p:', d % P == 0)
if d % P == 0:
    ripple(v, {33787: d // P})
    score(v, 'x_33787 set')

print('\n=== cone of x_30454 (drives x_12186) ===')
import collections
def bcone(roots_, maxn=200000):
    seen = set(roots_); q = collections.deque(roots_); free = set()
    while q:
        u = q.popleft()
        a = definer.get(u)
        if a is None:
            free.add(u); continue
        for w in set(x for m in polys[a] for x in m):
            if w != u and w not in seen:
                seen.add(w); q.append(w)
    return seen, free
s, f = bcone([30454])
print('  cone size', len(s), 'free inputs', len(f))
print('  x_12186 =', v0[12186], '\n  x_30454 =', v0[30454], '\n  x_33612(free) =', v0[33612])
print('  x_12186 == x_30454 + p*x_33612 :', v0[12186] == v0[30454] + P * v0[33612])
print('  need x_30454 = H1618 mod p ; delta mod p =', (v0[30454] - H1618) % P)
print('  atom 31928 src:', src[31928])
for u in (21608, 34575):
    print(f'   x_{u} def={definer.get(u)} src={src[definer[u]][:120] if u in definer else "FREE"}')
