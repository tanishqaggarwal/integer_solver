"""jm step 19: the A = 0 configuration, priced exactly.

The seven residual atoms are
  a22229 = x_7068 - x_2099 - 7376877*x_642      a22230 = x_28730 - p*x_9413
  a35758 = x_29854 - p*x_1329                   a35759 = 5113045*x_7075*x_9118 - x_29854
  a35760 = x_31864 - p*x_10903                  a35761 = x_7075*x_8731 + x_31864
  a35762 = x_642 - p*x_17325
so A = 0  <=>  p | x_9118,  p | x_8731,  x_28730 = p*x_9413,  x_7068 = x_2099 + 7376877*p*x_17325.
A = 0 satisfies all twelve equations, so the whole instance then scores
39033 - out12.  Set it directly (x_8731 = x_9118 = 0 is the cleanest choice, and
both are zero-collateral inputs) and measure the price of the two congruences.
"""
import os, sys, json, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
import jm_05_engine as EN
from jm_14_sel import build
P = J.P
W = J.base_state()
LOG = '/home/user/integer_solver/solve_lab/s10/jm_azero.jsonl'
BEST = '/home/user/integer_solver/solve_lab/s10/jm_best.json'


def full(v):
    av = L.all_atom_values(v)
    fq = set(L.failing_eqs(av))
    nz = [a for a in range(L.NA) if av[a]]
    return len(fq - J.E12), L.NEQ - len(fq), nz, sorted(fq & J.E12)


def report(v, tag, f=None):
    o, s, nz, ins = full(v)
    sev = [a for a in J.SEVEN if a in nz]
    print(f'  {tag:<46} out12={o:<4} in12={len(ins):<3} score={s:<6} '
          f'A!=0:{sev}  other={[a for a in nz if a not in J.SS][:10]}',
          flush=True)
    if f:
        f.write(json.dumps({'tag': tag, 'out12': o, 'in12': len(ins),
                            'score': s, 'nonzero': nz}) + '\n')
        f.flush()
    return o, s, nz


def azero(v, tracker=True):
    """force A = 0"""
    v[8731] = 0
    v[9118] = 0
    v[1329] = 0
    v[10903] = 0
    v[17325] = 0
    v[9413] = 0
    J.fwd2(v, 2)
    v[642] = 0
    v[29854] = 0
    v[31864] = 0
    v[28730] = 0
    v[7068] = v[2099]
    J.fwd2(v, 2)
    if tracker:
        for _ in range(3):
            for tr, td in ((14853, 1308), (24548, 25442), (14623, 27522)):
                d = v[td] - W[td]
                if d and (v[tr] - W[tr]) != d:
                    v[tr] = W[tr] + d
            J.fwd2(v, 2)
            v[7068] = v[2099]
            J.fwd2(v, 2)
    return v


if __name__ == '__main__':
    f = open(LOG, 'a')
    print('control:')
    report(W, 'witness', f)
    best = (39026, None, 'witness')
    for nm, v0, tk in (('witness', list(W), True),
                       ('witness-notrack', list(W), False),
                       ('sel(0,1)', build(0, 1), True),
                       ('sel(0,0)', build(0, 0), True),
                       ('sel(1,1)', build(1, 1), True)):
        print(f'\n=== A=0 from {nm} (tracker={tk}) ===')
        v = azero(list(v0), tk)
        o, s, nz = report(v, f'{nm}: A=0 raw', f)
        vr, o2, s2, nz2 = EN.repair(list(v), verbose=False, maxit=30)
        for _ in range(3):
            vr = azero(vr, tk)
            vr, o2, s2, nz2 = EN.repair(list(vr), verbose=False, maxit=20)
        o3, s3, nz3 = report(vr, f'{nm}: A=0 + repair', f)
        if s3 > best[0]:
            best = (s3, vr, nm)
    print(f'\nBEST {best[0]} from {best[2]}')
    if best[0] > 39026 and best[1] is not None:
        json.dump({f'x_{i}': best[1][i] for i in range(L.NVARS)
                   if best[1][i] != 0}, open(BEST, 'w'))
        print('saved jm_best.json')
    f.close()
