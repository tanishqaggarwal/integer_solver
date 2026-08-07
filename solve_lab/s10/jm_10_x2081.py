"""jm step 10: the x_2081 multiplier route to congruence 1.

a3576 = x_2081*(x_6418 - C) - 15804267*x_26777,   x_26777 = p*x_3387.
If x_2081 == 0 mod p the whole product dies mod p and x_3387 can absorb the rest
over Z, so x_6418 -- hence x_2099 and C0 -- becomes free mod p at NO a3576 cost.
Price x_2081 at its p-quantised values and build the joint move.
"""
import os, sys, time, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
import jm_05_engine as EN
P = J.P
W = J.base_state()
R0 = J.resid(W)
C = -[c for m, c in L.polys[3576].items() if m == (2081,)][0]
M = 15804267
LOG = '/home/user/integer_solver/solve_lab/s10/jm_x2081.jsonl'


def show(v, tag, f=None):
    c, s, nz, av = EN.state(v)
    R = J.resid(v)
    mv = ''.join(n for n, a, b in zip('ABC', R, R0) if a != b)
    print(f'  {tag:<46} out12={c:<4} score={s:<6} moves={mv:<4} broken={nz[:12]}',
          flush=True)
    if f:
        f.write(json.dumps({'tag': tag, 'out12': c, 'score': s, 'moves': mv,
                            'broken': nz}) + '\n')
        f.flush()
    return c, s, nz


if __name__ == '__main__':
    f = open(LOG, 'a')
    print(f'x_2081 = {W[2081]}, C bits {C.bit_length()}, M = {M}')
    print(f'a29090 (defines x_2099): {L.atom_src[L.definer[2099]]}')
    for u in (37158, 25297):
        d = J.definer.get(u)
        print(f'  x_{u}: free={d is None} val={str(W[u])[:30]} '
              f'src={L.atom_src[d][:120] if d is not None else ""}')
    print('\n--- price x_2081 at p-quantised / small values ---')
    for val in (0, 1, -1, 2, P, 2 * P, -P, P + 1, M, M * P):
        v = list(W); v[2081] = val; J.fwd2(v, 2)
        show(v, f'x_2081 = {"p" if val == P else str(val)[:22]}', f)

    print('\n--- joint: x_2081 = p, x_6418 = C + M*k, x_3387 = k ---')
    for k in (1, 2, 7, 1000003):
        v = list(W)
        v[2081] = P
        v[6418] = C + M * k
        v[3387] = k
        J.fwd2(v, 2)
        c, s, nz = show(v, f'x_2081=p, k={k}', f)

    print('\n--- joint without touching x_3387 (let fwd solve it) ---')
    for k in (1, 1000003):
        v = list(W)
        v[2081] = P
        v[6418] = C + M * k
        J.fwd2(v, 2)
        show(v, f'x_2081=p, x_6418=C+M*{k}, x_3387 auto', f)

    print('\n--- with repair engine, keeping congruence 1 relaxed ---')
    for val in (P, 0):
        for k in (1, 1000003):
            v = list(W)
            v[2081] = val
            v[6418] = C + M * k
            v[3387] = k
            J.fwd2(v, 2)
            if not EN.keep_c1(v):
                print(f'   x_2081={val} k={k}: C0 NOT moved, skip')
                continue
            vr, c2, s2, nz2 = EN.repair(list(v), keep=EN.keep_c1, verbose=False)
            print(f'   x_2081={"p" if val == P else val} k={k}: repaired '
                  f'out12={c2} score={s2} broken={nz2[:12]}', flush=True)
            f.write(json.dumps({'tag': f'rep x2081={val} k={k}', 'out12': c2,
                                'score': s2, 'broken': nz2}) + '\n')
            f.flush()
    f.close()
