"""WR step 1: detach x_26064 (a37694 becomes a check), drive the whole wire to w,
forward-evaluate in that frame, measure the score exactly."""
import os, sys, collections, json, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
ROOT, ROOTATOM = 26064, 37694

# ---- frame: detach the wire root ------------------------------------------
DETACH = {ROOT: ROOTATOM}
definer = {t: a for t, a in L.definer.items() if t not in DETACH}
atom_out = {a: o for a, o in L.atom_out.items() if a not in set(DETACH.values())}
ORDER = [t for t in ad.ORDER if t not in DETACH]
FREE = set(t for t in range(L.NVARS) if t not in definer)
CHECKS = [a for a in range(L.NA) if a not in atom_out]


def fwd(v, rounds=8):
    for _ in range(rounds):
        for u in ORDER:
            nv = T.solve_lin(definer[u], u, v)
            if nv is not None:
                v[u] = nv
    return v


def report(v, tag):
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    fail = L.failing_eqs(av)
    sc = L.NEQ - len(fail)
    nzc = [a for a in nz if a not in atom_out]
    nzg = [a for a in nz if a in atom_out]
    print(f'{tag}: score={sc} failing={len(fail)} nonzero_atoms={len(nz)} '
          f'(checks {len(nzc)}, broken-gates {len(nzg)})', flush=True)
    return av, nz, fail, sc


if __name__ == '__main__':
    base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
    WIRE = sorted(u for u in range(L.NVARS) if base[u] == P)
    print(f'wire members with value p in the deliverable: {len(WIRE)}')
    print(f'root x_{ROOT} in wire: {ROOT in WIRE}; a{ROOTATOM} in {len(L.atom2eq[ROOTATOM])} eqs')
    print(f'a{ROOTATOM} src: {L.atom_src[ROOTATOM][:120]}')
    print(f'a{ROOTATOM} eqs: {sorted(L.atom2eq[ROOTATOM])}')

    # sanity: frame keeps the deliverable on-manifold
    b2 = list(base); fwd(b2)
    report(b2, 'deliverable in detached frame')

    RESULTS = {}
    vals = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [1, -1, 2, 0, 3, -2, P]
    for W in vals:
        v = list(base)
        for u in WIRE:
            v[u] = W
        v[ROOT] = W
        fwd(v, rounds=8)
        # re-pin any wire member that drifted (should not, copies are exact)
        held = sum(1 for u in WIRE if v[u] == W)
        av, nz, fail, sc = report(v, f'w={W} (wire held {held}/{len(WIRE)})')
        nzc = [a for a in nz if a not in atom_out]
        print(f'   nonzero checks: {sorted(nzc)[:40]}')
        print(f'   failing eqs   : {sorted(fail)[:40]}')
        RESULTS[W] = sc
        T.save(v, os.path.join(HERE, f'wr_w{W}.json'.replace('-', 'm')))
    print()
    for W, sc in RESULTS.items():
        print(f'   w={W:<20} score {sc}')
