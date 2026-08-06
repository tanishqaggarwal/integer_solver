"""jm step 18: exhaustive TWO-further-move search from a relaxed state.

From the C1 state (out12 = 13, only a3576 broken) the complete set of free inputs
that can touch a3576's thirteen equations is 32.  All 496 pairs of them, over a
delta grid that includes each atom-zeroing solve, are measured exactly.  That is
an exhaustive 3-parameter search around the cheapest congruence-1 relaxation.

usage: python3 jm_18_pair2.py <state> START END
"""
import os, sys, json, time, itertools
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
import jm_05_engine as EN
from jm_17_comp import state, candidates
from jm_14_sel import relax
P = J.P
W = J.base_state()
LOG = '/home/user/integer_solver/solve_lab/s10/jm_pair2.jsonl'
BEST = '/home/user/integer_solver/solve_lab/s10/jm_best.json'


def done_set():
    s = set()
    if os.path.exists(LOG):
        for ln in open(LOG):
            try:
                r = json.loads(ln)
                s.add((r['state'], r['u'], r['w'], r['i']))
            except Exception:
                pass
    return s


if __name__ == '__main__':
    name = sys.argv[1]
    a, b = int(sys.argv[2]), int(sys.argv[3])
    v0 = state(name)
    F, feqs, atoms, (o0, s0, nz0) = candidates(v0)
    r0 = relax(v0)
    PAIRS = list(itertools.combinations(F, 2))
    print(f'{name}: out12={o0}, {len(F)} inputs, {len(PAIRS)} pairs; '
          f'batch [{a},{b})', flush=True)
    have = done_set()
    f = open(LOG, 'a')
    t0 = time.time()
    best = o0
    G = [1, 1000003, P]

    def solves(u, v):
        out = []
        for at in nz0:
            if u in L.avars[at]:
                t = T.solve_lin(at, u, v)
                if t is not None:
                    out.append(t - v[u])
        return out

    for idx in range(a, min(b, len(PAIRS))):
        u, w = PAIRS[idx]
        du = G + solves(u, v0)
        dw = G + solves(w, v0)
        combos = [(x, y) for x in du for y in dw]
        for i, (x, y) in enumerate(combos):
            if (name, u, w, i) in have:
                continue
            v = list(v0)
            v[u] += x
            v[w] += y
            J.fwd2(v, 2)
            if relax(v) != r0:
                continue
            o, s, nz, av = EN.state(v)
            if o < best or o <= o0 - 2:
                f.write(json.dumps({'state': name, 'u': u, 'w': w, 'i': i,
                                    'dx': str(x), 'dy': str(y), 'out12': o,
                                    'score': s, 'broken': nz}) + '\n')
                f.flush()
            if o < best:
                best = o
                print(f'  * x_{u}+{str(x)[:10]} & x_{w}+{str(y)[:10]}: '
                      f'out12={o} score={s} broken={nz[:10]}', flush=True)
                if o < 7:
                    json.dump({f'x_{k}': v[k] for k in range(L.NVARS)
                               if v[k] != 0}, open(BEST, 'w'))
        if (idx - a) % 40 == 0:
            f.write(json.dumps({'state': name, 'u': u, 'w': w, 'i': -1,
                                'out12': -1}) + '\n')
            f.flush()
            print(f'  ..{idx}/{len(PAIRS)} best={best} '
                  f'({time.time()-t0:.0f}s)', flush=True)
    f.close()
    print(f'batch done best={best} ({time.time()-t0:.0f}s)', flush=True)
