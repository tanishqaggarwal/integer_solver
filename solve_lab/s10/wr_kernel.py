"""WR step 9: measure the KERNEL deformations (0 identity rows broken) in the
modern frame, and the constrained deformation that pins the four residual
multipliers to 1."""
import os, sys, json, collections, random
from fractions import Fraction
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W
import wr_rows as R
P = ad.P
WIRE, widx, rows, RE = R.WIRE, R.widx, R.rows, R.RE
N = len(WIRE)
F = W.F_WIRE
MULT = [17499, 22665, 28961, 28599]

base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
b2 = list(base); F.fwd(b2)

KB = json.load(open(os.path.join(HERE, 'wirekernel.json')))
assert KB['wire'] == WIRE, (len(KB['wire']), len(WIRE))
BAS = KB['basis']
# verify against my rows
for b in BAS:
    for e in RE:
        assert sum(c * b[j] for j, c in rows[e].items()) == 0, e
print('kernel basis verified against my 219 identity rows')


def apply(d, tag, rounds=10, run_engine=False):
    v = list(b2)
    for j, u in enumerate(WIRE):
        v[u] = P + d[j]
    F.fwd(v, rounds=rounds)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    held = sum(1 for j, u in enumerate(WIRE) if v[u] == P + d[j])
    print(f'{tag}: score={L.NEQ-len(fail)} failing={len(fail)} nonzero={len(nz)} '
          f'wire held {held}/{N}', flush=True)
    ident_bad = [e for e in RE if e in set(fail)]
    print(f'   failing identity rows: {len(ident_bad)}; '
          f'multipliers now: {[(u, str(v[u])[:14]) for u in MULT]}')
    return v, nz, fail


if __name__ == '__main__':
    for i, b in enumerate(BAS):
        apply(b, f'kernel basis {i}')
    rnd = random.Random(3)
    for t in range(2):
        d = [sum(rnd.randrange(-3, 4) * b[j] for b in BAS) for j in range(N)]
        apply(d, f'random kernel combo {t}')
    # kernel + uniform
    for c in (1,):
        d = [c + BAS[0][j] for j in range(N)]
        apply(d, f'uniform({c}) + kernel0')
