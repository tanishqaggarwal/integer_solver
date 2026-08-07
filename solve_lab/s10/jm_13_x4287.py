"""jm step 13: the x_4287 selector route to congruence 1.

x_2099 = x_37158 + x_25297,  x_25297 = x_21279 * x_9118,  x_21279 = f(x_4287).
At the witness x_4287 = 0 so x_25297 = 0 and C0 = x_7068 - x_37158 is rigid.
At x_4287 = 1, x_21279 = 1 and C0 = x_7068 - x_37158 - x_9118, so x_9118 -- a
ZERO-collateral free input -- becomes a knob on C0 mod p.

Turning the selector on activates two more constant pins,
   a3568 = x_4287*(x_31861 - K5),  a3570 = x_4287*(x_14865 - K6),
whose pinned variables x_31861 and x_14865 are THEMSELVES free inputs, so they
can be set to the constants.  That is a 4-parameter joint move; no single-move
scan can see it.
"""
import os, sys, json, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
import jm_05_engine as EN
P = J.P
W = J.base_state()
R0 = J.resid(W)
LOG = '/home/user/integer_solver/solve_lab/s10/jm_x4287.jsonl'


def const_of(a, u):
    """atom a = x_sel * (x_u - K); return K"""
    for m, c in L.polys[a].items():
        if len(m) == 1:
            return -c
    return None


def show(v, tag, f=None):
    c, s, nz, av = EN.state(v)
    c1 = (v[7068] - v[2099]) % P != (W[7068] - W[2099]) % P
    c2 = (v[28730] - W[28730]) % P != 0
    print(f'  {tag:<50} out12={c:<4} score={s:<6} C1={int(c1)} C2={int(c2)} '
          f'broken={nz[:12]}', flush=True)
    if f:
        f.write(json.dumps({'tag': tag, 'out12': c, 'score': s, 'C1': c1,
                            'C2': c2, 'broken': nz}) + '\n')
        f.flush()
    return c, s, nz, c1, c2


if __name__ == '__main__':
    f = open(LOG, 'a')
    K5 = const_of(3568, 31861)
    K6 = const_of(3570, 14865)
    print(f'a3568 src: {L.atom_src[3568][:120]}')
    print(f'a3570 src: {L.atom_src[3570][:120]}')
    print(f'K5 bits {K5.bit_length()}   x_31861 = {W[31861]}')
    print(f'K6 bits {K6.bit_length()}   x_14865 = {W[14865]}')
    print(f'x_21279 = {W[21279]}   x_9118 = {W[9118]}  (== p: {W[9118] == P})')
    print(f'x_25297 = {W[25297]}   x_37158 = {str(W[37158])[:40]}')
    print()

    print('--- selector on, pins fed ---')
    v = list(W); v[4287] = 1; J.fwd2(v, 2)
    show(v, 'x_4287=1', f)
    v = list(W); v[4287] = 1; v[31861] = K5; J.fwd2(v, 2)
    show(v, 'x_4287=1, x_31861=K5', f)
    v = list(W); v[4287] = 1; v[14865] = K6; J.fwd2(v, 2)
    show(v, 'x_4287=1, x_14865=K6', f)
    v0 = list(W); v0[4287] = 1; v0[31861] = K5; v0[14865] = K6; J.fwd2(v0, 2)
    show(v0, 'x_4287=1, both pins fed', f)
    print(f'   x_21279 now {v0[21279]}  x_25297 {str(v0[25297])[:30]}')

    print('\n--- + repair ---')
    vr, c, s, nz = EN.repair(list(v0), verbose=True, tag='x4287', maxit=30)
    show(vr, 'x_4287=1 both pins + repair', f)

    print('\n--- then move x_9118 (C0 knob) ---')
    for base, bt in ((v0, 'pins-fed'), (vr, 'repaired')):
        for d in (1, 1000003):
            v = list(base); v[9118] = v[9118] + d; J.fwd2(v, 2)
            c, s, nz, c1, c2 = show(v, f'{bt} + x_9118+{d}', f)
            if c1:
                v2, c2_, s2_, nz2 = EN.repair(list(v), keep=EN.keep_c1,
                                              verbose=False, maxit=30)
                show(v2, f'{bt} + x_9118+{d} + repair', f)

    print('\n--- x_4287=1 also unpins?  price x_6418 / x_28730 from there ---')
    for u, d in ((6418, 1000003), (28730, 1000003), (7068, 1000003)):
        v = list(vr); v[u] = v[u] + d
        J.fwd2(v, 2)
        for tr, td in ((14853, 1308), (24548, 25442), (14623, 27522)):
            v[tr] += v[td] - vr[td]
        J.fwd2(v, 2)
        show(v, f'repaired-x4287 + x_{u}+{d} + trackers', f)
    f.close()
