"""jm step 14: all four selector settings with every constant pin fed, then a
broad single-input rescan of the leftovers (recording EVERY measurement, not
just improvements, so lateral 2-step paths are visible).

selectors:  x_2081 (pins a3576: x_6418 == C4,  a3578: x_12553 == C3)
            x_4287 (pins a3568: x_31861 == K5, a3570: x_14865 == K6)
CHUNKED: python3 jm_14_sel.py <phase>   phase in {sel, scan}
"""
import os, sys, json, time, itertools
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
import jm_05_engine as EN
P = J.P
W = J.base_state()
LOG = '/home/user/integer_solver/solve_lab/s10/jm_sel.jsonl'
STATE = '/home/user/integer_solver/solve_lab/s10/jm_sel_best.json'


def const_of(a):
    for m, c in L.polys[a].items():
        if len(m) == 1:
            return -c
    return None


C4 = const_of(3576)      # x_6418 pin
C3 = const_of(3578)      # x_12553 pin
K5 = const_of(3568)      # x_31861 pin
K6 = const_of(3570)      # x_14865 pin
PINVAL = {6418: C4, 12553: C3, 31861: K5, 14865: K6}


def relax(v):
    return ((v[7068] - v[2099]) % P != (W[7068] - W[2099]) % P,
            (v[28730] - W[28730]) % P != 0)


def show(v, tag, f=None):
    c, s, nz, av = EN.state(v)
    r1, r2 = relax(v)
    print(f'  {tag:<48} out12={c:<4} score={s:<6} C1={int(r1)} C2={int(r2)} '
          f'broken={nz[:12]}', flush=True)
    if f:
        f.write(json.dumps({'tag': tag, 'out12': c, 'score': s, 'C1': r1,
                            'C2': r2, 'broken': nz}) + '\n')
        f.flush()
    return c, s, nz, r1, r2


def close_frame(v, rounds=2):
    J.fwd2(v, rounds)
    for _ in range(3):
        for tr, td in ((14853, 1308), (24548, 25442), (14623, 27522)):
            d = v[td] - W[td]
            if d and (v[tr] - W[tr]) != d:
                v[tr] = W[tr] + d
        J.fwd2(v, rounds)
    for a, h, hin in ((29539, 29967, 30163), (7930, 7927, 11052),
                      (21617, 36864, 5040), (3576, 26777, 3387)):
        tgt = T.solve_lin(a, h, v)
        if tgt is None:
            continue
        d = J.definer.get(h)
        if d is None:
            continue
        vv = list(v); vv[h] = tgt
        nv = T.solve_lin(d, hin, vv)
        if nv is not None and nv != v[hin]:
            v[hin] = nv
            J.fwd2(v, rounds)
    return v


def build(s2081, s4287, extra=None):
    v = list(W)
    v[2081] = s2081
    v[4287] = s4287
    if s2081 == 1:
        v[6418] = C4
        v[12553] = C3
    if s4287 == 1:
        v[31861] = K5
        v[14865] = K6
    if extra:
        for u, d in extra.items():
            v[u] = v[u] + d
    close_frame(v)
    return v


if __name__ == '__main__':
    phase = sys.argv[1] if len(sys.argv) > 1 else 'sel'
    f = open(LOG, 'a')
    if phase == 'sel':
        best = None
        for s1 in (0, 1):
            for s2 in (0, 1):
                print(f'\n=== x_2081={s1}  x_4287={s2} ===')
                v = build(s1, s2)
                show(v, f'sel({s1},{s2}) pins fed', f)
                vr, c, s, nz = EN.repair(list(v), verbose=False, maxit=30)
                c, s, nz, r1, r2 = show(vr, f'sel({s1},{s2}) + repair', f)
                # add the congruence drivers on top
                for tag, ex in (('+x_9118', {9118: 1000003}),
                                ('+x_6418', {6418: 1000003}),
                                ('+x_28730', {28730: 1000003}),
                                ('+x_9118+x_28730', {9118: 1000003, 28730: 1000003}),
                                ('+x_6418+x_28730', {6418: 1000003, 28730: 1000003})):
                    v2 = list(vr)
                    for u, d in ex.items():
                        v2[u] = v2[u] + d
                    close_frame(v2)
                    c2, s2_, nz2, q1, q2 = show(v2, f'sel({s1},{s2}) {tag}', f)
                    if q1 or q2:
                        v3, c3, s3, nz3 = EN.repair(
                            list(v2), keep=lambda vv: relax(vv) == (q1, q2),
                            verbose=False, maxit=30)
                        c3, s3, nz3, q1b, q2b = show(
                            v3, f'sel({s1},{s2}) {tag} + repair', f)
                        if q1b and q2b and (best is None or c3 < best[0]):
                            best = (c3, s3, v3, f'sel({s1},{s2}) {tag}')
                        if q1b and q2b and c3 < 7:
                            json.dump({f'x_{i}': v3[i] for i in range(L.NVARS)
                                       if v3[i] != 0}, open(STATE, 'w'))
        if best:
            print(f'\nBEST both-congruence move: out12={best[0]} score={best[1]} '
                  f'via {best[3]}')
            json.dump({f'x_{i}': best[2][i] for i in range(L.NVARS)
                       if best[2][i] != 0}, open(STATE, 'w'))
    f.close()
