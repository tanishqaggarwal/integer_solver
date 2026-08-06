"""Step 20: (a) mod-p cancellation test at stateC1; (b) alternative defect placements."""
import sys, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9/quad')
from common import *

CODES, _ = H.load_equations()
H688 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
H1618 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
dd = pickle.load(open('atoms.pkl', 'rb')); eq_terms = dd['eq_terms']


def score(v, tag, quiet=False):
    nz = nz_checks(v); ng = nz_gates(v)
    live = [a for a in range(len(polys)) if evalpoly(polys[a], v) != 0]
    ff = H.evaluate(CODES, v, eqs_of(live))
    if not quiet:
        print(f'[{tag:36s}] nz={len(nz):3d} ng={len(ng)} FAIL={len(ff):4d} {sorted(set(nz+ng))}')
    return nz, ng, ff


v = H.load_assignment('quad/stateC1.json')
nz, ng, ff = score(v, 'stateC1')
A = evalpoly(polys[31670], v); B = evalpoly(polys[31672], v)
print(f'A(31670) = {A}\nB(31672) = {B}')
print('A mod p =', A % P, '\nB mod p =', B % P)

# handle structure -> what shifts are available
for u in (29309, 36358, 3915, 105):
    d = definer.get(u)
    print(f'  x_{u}: {"FREE" if d is None else src[d][:90]}  val={v[u]}')
print('  x_3915 == p:', v[3915] == P)

print('\n--- per failing equation: value = m*(c1*A + c2*B) ---')
for i in ff:
    m, sq, tl = eq_terms[i]
    c1 = sum(c for c, a in tl if a == 31670)
    c2 = sum(c for c, a in tl if a == 31672)
    val = H.resid(CODES, v, i)
    ok = (c1 * (A % P) + c2 * (B % P)) % P == 0
    print(f'  eq {i:6d}: m={m:3d} sq={int(sq)} c1={c1:4d} c2={c2:4d}  '
          f'c1*A+c2*B = 0 mod p ? {ok}   |val|={abs(val).bit_length()}b')

print('\n=== ALTERNATIVE PLACEMENTS ===')
vA2 = H.load_assignment('quad/stateA2.json')
PIN33462 = 97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
PIN22152 = 126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506

# from stateC1, revert x_33462 to its pinned residue (keeps 31672) -> chain-1 defect moves to 33929
v2 = H.load_assignment('quad/stateC1.json')
ripple(v2, {33462: PIN33462})
score(v2, 'C1, x_33462 restored (31672 ok)')

# from stateC1, revert x_22152 to pinned (keeps 31670) -> chain-2 defect moves to 2423
v3 = H.load_assignment('quad/stateC1.json')
ripple(v3, {22152: PIN22152})
score(v3, 'C1, x_22152 restored (31670 ok)')

# both restored
v4 = H.load_assignment('quad/stateC1.json')
ripple(v4, {33462: PIN33462, 22152: PIN22152})
score(v4, 'C1, both restored')
