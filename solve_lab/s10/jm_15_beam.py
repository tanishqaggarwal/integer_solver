"""jm step 15: beam search over repair moves from a relaxed state.

The greedy engine plateaus because closing one pin opens another (the residual is
conserved).  A beam keeps the K best states by (out12, #broken) and expands each
by every depth-3 repair candidate for every broken atom, plus explicit TRACKER
moves for every pin-shaped atom  K*(free - derived) - handle.

CHUNKED: python3 jm_15_beam.py <state> <levels> [beam]
state in {C1_6418, C2_28730, SEL01, SEL01_C2, SEL00, SEL11}
Every visited state is appended to jm_beam.jsonl.
"""
import os, sys, json, time, itertools, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
import jm_05_engine as EN
from jm_14_sel import build, close_frame, relax, C3, C4, K5, K6
P = J.P
W = J.base_state()
LOG = '/home/user/integer_solver/solve_lab/s10/jm_beam.jsonl'
BEST = '/home/user/integer_solver/solve_lab/s10/jm_best.json'


def start(name):
    if name == 'C1_6418':
        v = list(W); v[6418] += 1000003; J.fwd2(v, 2)
    elif name == 'C2_28730':
        v = list(W); v[28730] += 1000003; J.fwd2(v, 2)
        v[24548] += v[25442] - W[25442]; J.fwd2(v, 2)
    elif name == 'SEL01':
        v = build(0, 1)
    elif name == 'SEL01_C2':
        v = build(0, 1, {28730: 1000003})
    elif name == 'SEL00':
        v = build(0, 0)
    elif name == 'SEL11':
        v = build(1, 1)
    else:
        raise SystemExit('unknown state')
    return v


def key(v):
    c, s, nz, av = EN.state(v)
    return (c, len(nz)), c, s, nz


def cands(v, nz, depth=3):
    out = []
    for a in nz:
        out += EN.gen(a, v, depth)
        # tracker moves: for every var of a that is free, set it so a vanishes
        for u in sorted(L.avars[a]):
            if u not in J.FREESET:
                continue
            t = T.solve_lin(a, u, v)
            if t is not None and t != v[u]:
                out.append((u, t))
    seen = set()
    res = []
    for u, nv in out:
        if (u, nv) in seen:
            continue
        seen.add((u, nv))
        res.append((u, nv))
    return res


if __name__ == '__main__':
    name = sys.argv[1]
    levels = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    BW = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    v0 = start(name)
    k0, c0, s0, nz0 = key(v0)
    r0 = relax(v0)
    print(f'{name}: start out12={c0} score={s0} C=({int(r0[0])},{int(r0[1])}) '
          f'broken={nz0}', flush=True)
    f = open(LOG, 'a')
    beam = [(k0, v0, nz0, ())]
    seen = {tuple(sorted(nz0))}
    best = (c0, s0, v0, ())
    t0 = time.time()
    for lv in range(levels):
        nxt = []
        for k, v, nz, path in beam:
            for u, nv in cands(v, nz):
                tr = list(v)
                tr[u] = nv
                J.fwd2(tr, 2)
                r = relax(tr)
                if r != r0:                      # must keep the relaxation
                    continue
                kk, cc, ss, nn = key(tr)
                sig = tuple(sorted(nn))
                if sig in seen:
                    continue
                seen.add(sig)
                nxt.append((kk, tr, nn, path + (u,)))
                f.write(json.dumps({'state': name, 'level': lv, 'path': list(path + (u,)),
                                    'out12': cc, 'score': ss, 'broken': nn}) + '\n')
                if cc < best[0]:
                    best = (cc, ss, tr, path + (u,))
                    print(f'  L{lv} NEW BEST out12={cc} score={ss} via '
                          f'{path + (u,)} broken={nn}', flush=True)
        f.flush()
        nxt.sort(key=lambda t: t[0])
        beam = nxt[:BW]
        print(f'  level {lv}: {len(nxt)} states, beam now '
              f'{[t[0] for t in beam]} ({time.time()-t0:.0f}s)', flush=True)
        if not beam:
            break
    f.close()
    print(f'{name} BEST out12={best[0]} score={best[1]} path={best[3]}')
    if best[0] < 7 and all(relax(best[2])):
        json.dump({f'x_{i}': best[2][i] for i in range(L.NVARS)
                   if best[2][i] != 0}, open(BEST, 'w'))
        print('saved jm_best.json')
