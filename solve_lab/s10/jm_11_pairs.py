"""jm step 11: exact-cancellation PAIR sweep over the residue-moving free inputs.

For a pair (u, w) and a functional phi in {RA, RB, RC, A1} that both move, set
    x_u += delta,   x_w += -delta * phi(u) / phi(w)   (mod p)
so phi is left fixed while the other functionals generically move.  That is the
1-codimension slice of the affine subspace; every point of it is measured
EXACTLY (fwd2 + failing_eqs), never predicted.

CHUNKED + RESUMABLE:  python3 jm_11_pairs.py START END
Records append to jm_pairs.jsonl keyed by (u,w,phi,delta).
"""
import os, sys, json, time, itertools, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
import jm_05_engine as EN
from jm_07_lin import load
P = J.P
LOG = '/home/user/integer_solver/solve_lab/s10/jm_pairs.jsonl'

MOVERS = [91, 438, 490, 1203, 1530, 1603, 2081, 2441, 2455, 2460, 3420, 3545,
          3591, 4287, 4701, 5090, 5643, 5910, 6312, 6418, 6821, 7142, 7401,
          8370, 8487, 9260, 10513, 10836, 11189, 11368, 12054, 12095, 12553,
          13195, 13710, 14808, 14823, 14853, 16348, 16586, 16827, 17083, 17378,
          17406, 17760, 18022, 20564, 20837, 21074, 21103, 21112, 21702, 22351,
          22513, 22562, 22623, 23262, 23751, 24365, 24517, 24548, 26005, 27393,
          28005, 28645, 28713, 29785, 30448, 31899, 33287, 34175, 34974, 35179,
          35740, 35839, 35979, 36044, 37147, 37862, 38055, 38625, 7068, 28730]


def phis():
    cols, dR, dA1 = load()
    out = {}
    for u in MOVERS:
        if u not in dR:
            continue
        out[u] = {'RA': dR[u][0], 'RB': dR[u][1], 'RC': dR[u][2],
                  'A1': dA1[u]}
    return out


def jobs(PH):
    js = []
    for u, w in itertools.combinations(sorted(PH), 2):
        for ph in ('RA', 'RB', 'RC'):
            if PH[u][ph] and PH[w][ph]:
                js.append((u, w, ph))
    return js


def done_set():
    s = set()
    if os.path.exists(LOG):
        for ln in open(LOG):
            try:
                r = json.loads(ln)
                s.add((r['u'], r['w'], r['phi'], r['delta']))
            except Exception:
                pass
    return s


if __name__ == '__main__':
    a, b = int(sys.argv[1]), int(sys.argv[2])
    PH = phis()
    JS = jobs(PH)
    print(f'{len(JS)} (pair, functional) jobs; batch [{a},{b})', flush=True)
    W = J.base_state()
    R0 = J.resid(W)
    have = done_set()
    f = open(LOG, 'a')
    t0 = time.time()
    DELTAS = [1, 1000003]
    for i in range(a, min(b, len(JS))):
        u, w, ph = JS[i]
        r = -PH[u][ph] * pow(PH[w][ph], -1, P) % P
        for d in DELTAS:
            key = (u, w, ph, d)
            if key in have:
                continue
            e = d * r % P
            v = list(W)
            v[u] += d
            v[w] += e
            J.fwd2(v, 2)
            c, s, nz, av = EN.state(v)
            R = J.resid(v)
            mv = ''.join(n for n, x, y in zip('ABC', R, R0) if x != y)
            a1 = (v[28730] - W[28730]) % P != 0
            f.write(json.dumps({'u': u, 'w': w, 'phi': ph, 'delta': d,
                                'ratio': str(r), 'out12': c, 'score': s,
                                'moves': mv, 'A1': a1, 'broken': nz}) + '\n')
            f.flush()
            if c <= 16:
                print(f'  x_{u}+{d} & x_{w} (kill {ph}): out12={c} score={s} '
                      f'moves={mv} A1={a1} broken={nz[:10]}', flush=True)
        if (i - a) % 50 == 0:
            print(f'  ..{i} ({time.time()-t0:.0f}s)', flush=True)
    f.close()
    print(f'batch done ({time.time()-t0:.0f}s)', flush=True)
