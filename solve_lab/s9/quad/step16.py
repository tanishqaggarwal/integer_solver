"""Step 16: the terminal pins 31670 / 31672 and the double-knob repair."""
import sys, collections
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
    print(f'[{tag:36s}] nz={len(nz):3d} ng={len(ng)} FAIL={len(ff):4d} {sorted(set(nz+ng))[:20]}')
    return nz, ng, ff


for a in (31670, 31672, 33929, 2423, 26731):
    print(f'atom {a}: {src[a][:230]}')
print()
for u in (29309, 29311, 33462, 22152, 8778, 33787, 32253, 14768):
    d = definer.get(u)
    print(f'x_{u}: {"FREE" if d is None else src[d][:110]}  val_bits={v0[u].bit_length()}')

# handle chains
def pchain(u, n=6):
    out = []
    for _ in range(n):
        d = definer.get(u)
        if d is None:
            out.append(f'x_{u}=FREE'); break
        out.append(f'x_{u}<-{src[d][:60]}')
        vs = [w for w in set(x for m in polys[d] for x in m) if w != u]
        if len(vs) != 1:
            out.append(f'   vars {vs}'); break
        u = vs[0]
    return ' | '.join(out)

print('\n29309 chain:', pchain(29309))
print('29311 chain:', pchain(29311))

print('\n=== double-knob test for the 688 chain ===')
k = (H688 - v0[19083]) % P
v = list(v0)
ripple(v, {8778: v0[8778] + k, 33462: v0[33462] + k})
print('x_19083 == H688 mod p:', (v[19083] - H688) % P == 0)
print('atom 33929 resid:', evalpoly(polys[33929], v) != 0, ' 31672 resid:', evalpoly(polys[31672], v) != 0)
score(v, 'x_8778,x_33462 += k')
# close 26731 with x_16742 and the p-handles
ripple(v, {16742: v[19083]})
score(v, ' + x_16742 := x_19083')
d = 8863713 * (v[18956] - H688)
print('688 numerator / p exact:', d % P == 0)
if d % P == 0:
    ripple(v, {7497: d // P})
    score(v, ' + x_7497 handle')
d = 6348691 * (v[8778] - v[16144])
print('33929 numerator/p exact:', d % P == 0)
if d % P == 0:
    ripple(v, {32253: d // P}); score(v, ' + x_32253 handle')
H.save_assignment(v, 'quad/stateB1.json')
