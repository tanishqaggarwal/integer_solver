"""jm step 16: the decisive experiment -- price the ALL-ZERO atom vector.

If every one of the seven residual atoms is zero, all twelve equations hold and
the score is 39033 - out12.  So the whole question "is the 7-equation prize
reachable" is exactly "can the seven be zeroed with out12 < 7".
Zero them through their detached variables (frame 2) from several starting
states, repair, and measure.
"""
import os, sys, json, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
import jm_05_engine as EN
from jm_14_sel import build, close_frame, relax
P = J.P
W = J.base_state()
LOG = '/home/user/integer_solver/solve_lab/s10/jm_zero.jsonl'
BEST = '/home/user/integer_solver/solve_lab/s10/jm_best.json'
PAIR = [(22229, 7068), (22230, 28730), (35758, 29854), (35761, 31864),
        (35762, 642)]


def zero_seven(v, rounds=4):
    for _ in range(rounds):
        for a, u in PAIR:
            t = T.solve_lin(a, u, v)
            if t is not None and t != v[u]:
                v[u] = t
        J.fwd2(v, 2)
    return v


def full(v):
    av = L.all_atom_values(v)
    f = set(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if av[a]]
    return len(f - J.E12), L.NEQ - len(f), nz, sorted(f & J.E12)


def report(v, tag, f=None):
    o, s, nz, ins = full(v)
    print(f'  {tag:<44} out12={o:<4} in12={len(ins):<3} score={s:<6} '
          f'nonzero={nz[:12]}', flush=True)
    if f:
        f.write(json.dumps({'tag': tag, 'out12': o, 'in12': len(ins),
                            'score': s, 'nonzero': nz}) + '\n')
        f.flush()
    return o, s, nz


if __name__ == '__main__':
    f = open(LOG, 'a')
    print('base:')
    report(W, 'witness', f)
    STARTS = [('witness', list(W)),
              ('sel(0,0)', build(0, 0)),
              ('sel(0,1)', build(0, 1)),
              ('sel(1,1)', build(1, 1)),
              ('sel(1,0)+x28730', build(1, 0, {28730: 1000003}))]
    best = None
    for nm, v0 in STARTS:
        print(f'\n=== start {nm} ===')
        report(v0, f'{nm} before', f)
        v = zero_seven(list(v0))
        o, s, nz = report(v, f'{nm} + zero seven', f)
        vr, o2, s2, nz2 = EN.repair(list(v), verbose=False, maxit=30)
        # re-zero after repair, then repair again
        for _ in range(3):
            vr = zero_seven(vr)
            vr, o2, s2, nz2 = EN.repair(list(vr), verbose=False, maxit=20)
        o3, s3, nz3 = report(vr, f'{nm} + zero + repair', f)
        if best is None or s3 > best[0]:
            best = (s3, vr, nm)
    print(f'\nBEST overall score {best[0]} from {best[2]}')
    if best[0] > 39026:
        json.dump({f'x_{i}': best[1][i] for i in range(L.NVARS)
                   if best[1][i] != 0}, open(BEST, 'w'))
        print('saved jm_best.json')
    f.close()
