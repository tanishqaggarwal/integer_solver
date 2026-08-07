"""WR step 10: frame with the WHOLE WIRE detached (every wire member is a free
parameter, every wire copy atom is a check).  Now a kernel deformation is
actually realisable.  Measure what it costs in the real circuit."""
import os, sys, json, collections, random
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
MULT = [17499, 22665, 28961, 28599]

DET = dict(W.F3)
for u in WIRE:
    if u in L.definer:
        DET[u] = L.definer[u]
FW = W.Frame(DET)
print(f'detached {len(DET)} variables (wire members with a definer: '
      f'{sum(1 for u in WIRE if u in L.definer)})')
print(f'free params {len(FW.FREE)}, checks {len(FW.CHECKS)}')

base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
b2 = list(base); FW.fwd(b2)
FW.report(b2, 'deliverable in wire-free frame')


def apply(d, tag, rounds=10):
    v = list(b2)
    for j, u in enumerate(WIRE):
        v[u] = P + d[j]
    FW.fwd(v, rounds=rounds)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    held = sum(1 for j, u in enumerate(WIRE) if v[u] == P + d[j])
    ident_bad = sorted(set(fail) & set(RE))
    print(f'{tag}: score={L.NEQ-len(fail)} failing={len(fail)} nonzero={len(nz)} '
          f'held {held}/{N}  identity-rows-failing {len(ident_bad)}', flush=True)
    return v, nz, fail


if __name__ == '__main__':
    KB = json.load(open(os.path.join(HERE, 'wirekernel.json')))
    BAS = KB['basis']
    apply([0] * N, 'd = 0 (control)')
    apply([1] * N, 'uniform d = 1')
    for i, b in enumerate(BAS):
        v, nz, fail = apply(b, f'kernel basis {i}')
        if i == 0:
            T.save(v, os.path.join(HERE, 'wr_kern0.json'))
            print(f'   nonzero atoms: {sorted(nz)[:60]}')
    rnd = random.Random(11)
    for t in range(3):
        d = [sum(rnd.randrange(-2, 3) * b[j] for b in BAS) for j in range(N)]
        apply(d, f'random kernel combo {t}')
