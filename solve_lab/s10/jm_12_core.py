"""jm step 12: focused multi-driver search in the congruence core.

DRIVERS are the only free inputs that can change C0 or A1 mod p:
    C0 = x_7068 - x_2099 :  x_6418, x_7068, x_2081, x_4287
    A1 = x_28730         :  x_28730
plus the inputs that reach the pins that then break (x_12553, x_13195).

For every driver subset of size <= 3 and every delta pattern we
  1. apply the drivers,
  2. CLOSE the frame: set the trackers x_14853 / x_24548 / x_14623 to track
     x_1308 / x_25442 / x_27522 exactly, and back-solve the p-handles
     x_30163 / x_11052 / x_5040 / x_3387,
  3. measure exactly, then run the repair engine,
and record which congruences the result actually relaxes.

CHUNKED + RESUMABLE:  python3 jm_12_core.py START END
"""
import os, sys, json, time, itertools, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
import jm_05_engine as EN
from jm_07_lin import load
P = J.P
LOG = '/home/user/integer_solver/solve_lab/s10/jm_core.jsonl'

DRIVERS = [6418, 7068, 2081, 4287, 28730, 12553, 13195, 14853, 24548, 14623,
           9118, 31861, 8731, 8976, 6467, 11099, 14865, 27711]
TRACK = [(14853, 1308), (24548, 25442), (14623, 27522)]
HANDLES = [(29539, 29967, 30163), (7930, 7927, 11052), (21617, 36864, 5040),
           (3576, 26777, 3387), (33796, 25295, 31339), (688, 14257, None)]
BOOL = {2081: [0, 1], 4287: [0, 1]}


def close_frame(v, W, rounds=2):
    """track the pins and back-solve every p-handle, then fwd."""
    J.fwd2(v, rounds)
    for _ in range(3):
        for tr, td in TRACK:
            d = v[td] - W[td]
            if d and (v[tr] - W[tr]) != d:
                v[tr] = W[tr] + d
        J.fwd2(v, rounds)
    # p-handles: solve the pin atom for the handle var, then back to its input
    for a, h, hin in HANDLES:
        if hin is None:
            continue
        tgt = T.solve_lin(a, h, v)
        if tgt is None:
            continue
        d = J.definer.get(h)
        if d is None:
            continue
        vv = list(v)
        vv[h] = tgt
        nv = T.solve_lin(d, hin, vv)
        if nv is not None and nv != v[hin]:
            v[hin] = nv
            J.fwd2(v, rounds)
    return v


def relaxations(v, W):
    c1 = (v[7068] - v[2099]) % P != (W[7068] - W[2099]) % P
    c2 = (v[28730] - W[28730]) % P != 0
    return c1, c2


def patterns():
    """(dict var->newvalue-offset-or-absolute, tag)"""
    out = []
    DS = [1, 1000003]
    for u in DRIVERS:
        for d in DS:
            out.append(({u: ('+', d)}, f'{u}+{d}'))
        if u in BOOL:
            for b in BOOL[u]:
                out.append(({u: ('=', b)}, f'{u}={b}'))
    for u, w in itertools.combinations(DRIVERS, 2):
        for d in DS:
            out.append(({u: ('+', d), w: ('+', d)}, f'{u}+{w}+{d}'))
            out.append(({u: ('+', d), w: ('+', P - d)}, f'{u}-{w}~{d}'))
        for u2, bs in BOOL.items():
            if u2 in (u, w):
                continue
        if u in BOOL:
            for b in BOOL[u]:
                for d in DS:
                    out.append(({u: ('=', b), w: ('+', d)}, f'{u}={b},{w}+{d}'))
    for u, w, z in itertools.combinations(DRIVERS, 3):
        out.append(({u: ('+', 1), w: ('+', 1), z: ('+', 1)}, f'{u},{w},{z}+1'))
        out.append(({u: ('+', 1000003), w: ('+', 1000003), z: ('+', 1000003)},
                    f'{u},{w},{z}+1000003'))
    return out


def done_set():
    s = set()
    if os.path.exists(LOG):
        for ln in open(LOG):
            try:
                s.add(json.loads(ln)['tag'])
            except Exception:
                pass
    return s


if __name__ == '__main__':
    a, b = int(sys.argv[1]), int(sys.argv[2])
    PT = patterns()
    print(f'{len(PT)} patterns; batch [{a},{b})', flush=True)
    W = J.base_state()
    have = done_set()
    f = open(LOG, 'a')
    t0 = time.time()
    best = 99
    for i in range(a, min(b, len(PT))):
        mv, tag = PT[i]
        if tag in have:
            continue
        v = list(W)
        for u, (op, d) in mv.items():
            v[u] = (v[u] + d) if op == '+' else d
        close_frame(v, W)
        c1, s1, nz1, av1 = EN.state(v)
        r1, r2 = relaxations(v, W)
        rec = {'tag': tag, 'raw': c1, 'raw_score': s1, 'c1': r1, 'c2': r2,
               'broken': nz1}
        if (r1 or r2) and c1 <= 30:
            keep = (lambda vv: relaxations(vv, W) == (r1, r2))
            vr, c2_, s2_, nz2 = EN.repair(list(v), keep=keep, verbose=False,
                                          maxit=25)
            rec.update({'rep': c2_, 'rep_score': s2_, 'rep_broken': nz2})
            if c2_ < best:
                best = c2_
                print(f'  * {tag}: raw={c1} rep={c2_} c1={r1} c2={r2} '
                      f'broken={nz2[:10]}', flush=True)
        f.write(json.dumps(rec) + '\n')
        f.flush()
        if (i - a) % 100 == 0:
            print(f'  ..{i} ({time.time()-t0:.0f}s)', flush=True)
    f.close()
    print(f'batch done ({time.time()-t0:.0f}s)', flush=True)
